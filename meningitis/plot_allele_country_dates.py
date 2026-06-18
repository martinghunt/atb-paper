#!/usr/bin/env python3
"""Plot country/year distribution for PASS meningococcal allele matches."""

from __future__ import annotations

import csv
import os
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


INPUTS = [
    ROOT / "WP_002226625.1.lexicmap.tsv.gz.perfect_match_samples.atb_meta.tsv",
    ROOT / "WP_002248791.1.lexicmap.tsv.gz.perfect_match_samples.atb_meta.tsv",
]

PLOT_PREFIX = ROOT / "supplementary_allele_country_dates"
PLOT_DATA = ROOT / "supplementary_allele_country_year_counts.tsv"
SUMMARY_TABLE = ROOT / "supplementary_allele_country_date_summary.tsv"
COUNTRY_SUMMARY_TABLE = ROOT / "supplementary_allele_country_summary.tsv"

MENINGOCOCCUS = "Neisseria meningitidis"


def allele_name(path: Path) -> str:
    return path.name.split(".lexicmap", 1)[0]


def normalize_country(value: str) -> str | None:
    country = (value or "").strip()
    if not country or country == "Not applicable":
        return None
    return country.split(":", 1)[0].strip()


def collection_year(value: str) -> int | None:
    years = [int(match) for match in re.findall(r"(?:19|20)\d{2}", value or "")]
    if not years:
        return None
    # Ranges such as 2006/2012 are plotted at the earliest plausible year,
    # which is the conservative choice for "present by year" statements.
    return min(years)


def date_format(value: str) -> str:
    date = (value or "").strip()
    if re.fullmatch(r"\d{4}", date):
        return "YYYY"
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", date):
        return "YYYY-MM-DD"
    if re.fullmatch(r"\d{4}/\d{4}", date):
        return "YYYY/YYYY"
    if date:
        return "partial_or_other"
    return "missing"


def load_records() -> tuple[pd.DataFrame, pd.DataFrame]:
    records: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []

    for path in INPUTS:
        allele = allele_name(path)
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

        pass_rows = [row for row in rows if row["hq_filter"] == "PASS"]
        nm_rows = [row for row in pass_rows if row["sylph_species"] == MENINGOCOCCUS]
        hq_counter = Counter(row["hq_filter"] for row in rows)
        species_counter = Counter(row["sylph_species"] for row in pass_rows)

        usable = []
        for row in nm_rows:
            country = normalize_country(row["country"])
            year = collection_year(row["collection_date"])
            if country is None or year is None:
                continue
            records.append(
                {
                    "allele": allele,
                    "sample_accession": row["sample_accession"],
                    "country_raw": row["country"],
                    "country": country,
                    "collection_date": row["collection_date"],
                    "collection_year": year,
                    "date_format": date_format(row["collection_date"]),
                }
            )
            usable.append((country, year))

        countries = sorted({country for country, _year in usable})
        years = [year for _country, year in usable]
        summaries.append(
            {
                "allele": allele,
                "total_rows": len(rows),
                "pass_rows": len(pass_rows),
                "pass_meningococcal_rows": len(nm_rows),
                "pass_non_meningococcal_rows": len(pass_rows) - len(nm_rows),
                "plotted_country_date_rows": len(usable),
                "plotted_countries": len(countries),
                "earliest_plotted_year": min(years) if years else "",
                "latest_plotted_year": max(years) if years else "",
                "hq_filter_counts": ";".join(
                    f"{key}:{value}" for key, value in sorted(hq_counter.items())
                ),
                "pass_species_counts": ";".join(
                    f"{key}:{value}" for key, value in sorted(species_counter.items())
                ),
            }
        )

    return pd.DataFrame(records), pd.DataFrame(summaries)


def write_tables(records: pd.DataFrame, summaries: pd.DataFrame) -> pd.DataFrame:
    if records.empty:
        raise SystemExit("No PASS meningococcal records with country and date were found.")

    grouped = (
        records.groupby(["allele", "country", "collection_year"], as_index=False)
        .agg(samples=("sample_accession", "nunique"))
        .sort_values(["allele", "country", "collection_year"])
    )
    grouped.to_csv(PLOT_DATA, sep="\t", index=False)
    country_summary = (
        records.groupby(["allele", "country"], as_index=False)
        .agg(
            samples=("sample_accession", "nunique"),
            earliest_year=("collection_year", "min"),
            latest_year=("collection_year", "max"),
        )
        .sort_values(["allele", "earliest_year", "country"])
    )
    country_summary.to_csv(COUNTRY_SUMMARY_TABLE, sep="\t", index=False)
    summaries.to_csv(SUMMARY_TABLE, sep="\t", index=False)
    return grouped


def make_plot(grouped: pd.DataFrame, summaries: pd.DataFrame) -> None:
    earliest_by_country = (
        grouped.groupby("country")["collection_year"].min().sort_values(kind="mergesort")
    )
    country_order = sorted(
        earliest_by_country.index,
        key=lambda country: (earliest_by_country[country], country),
    )
    y_positions = {country: index for index, country in enumerate(country_order)}
    alleles = list(summaries["allele"])
    max_count = int(grouped["samples"].max())
    marker_scale = 34

    fig_height = max(6.0, 0.34 * len(country_order) + 2.1)
    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(alleles),
        figsize=(11.5, fig_height),
        sharex=True,
        sharey=True,
    )
    fig.subplots_adjust(left=0.15, right=0.84, top=0.98, bottom=0.11, wspace=0.04)
    if len(alleles) == 1:
        axes = [axes]

    x_min = int(grouped["collection_year"].min()) - 2
    x_max = int(grouped["collection_year"].max()) + 1

    for ax, allele in zip(axes, alleles, strict=True):
        data = grouped[grouped["allele"] == allele].copy()
        sizes = 18 + marker_scale * (data["samples"] ** 0.5)
        ax.scatter(
            data["collection_year"],
            data["country"].map(y_positions),
            s=sizes,
            color="#4D4D4D",
            edgecolor="white",
            linewidth=0.7,
            alpha=0.88,
        )
        ax.grid(axis="x", color="#D9D9D9", linewidth=0.7)
        ax.grid(axis="y", color="#ECECEC", linewidth=0.5)
        ax.set_axisbelow(True)
        ax.set_xlim(x_min, x_max)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].set_yticks(range(len(country_order)), country_order)
    axes[0].invert_yaxis()
    for ax in axes:
        ax.tick_params(axis="y", length=0)

    axes[0].set_ylabel("Country")
    fig.supxlabel("Collection year")

    legend_counts = [count for count in [1, 5, 20, 100] if count <= max_count]
    legend_counts = sorted(set(legend_counts))
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#737373",
            markeredgecolor="white",
            markersize=((18 + marker_scale * (count**0.5)) ** 0.5),
            label=str(count),
        )
        for count in legend_counts
    ]
    fig.legend(
        handles=handles,
        title="Genomes",
        loc="center right",
        bbox_to_anchor=(0.98, 0.5),
        frameon=False,
        borderpad=0.2,
        labelspacing=0.9,
    )

    fig.savefig(f"{PLOT_PREFIX}.pdf", bbox_inches="tight")
    fig.savefig(f"{PLOT_PREFIX}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    records, summaries = load_records()
    grouped = write_tables(records, summaries)
    make_plot(grouped, summaries)
    print("Per-allele plotted counts:")
    for summary in summaries.itertuples(index=False):
        print(
            f"{summary.allele}: {summary.plotted_country_date_rows} genomes, "
            f"{summary.plotted_countries} countries, "
            f"{summary.earliest_plotted_year}-{summary.latest_plotted_year}"
        )
    print(f"Wrote {PLOT_PREFIX}.pdf")
    print(f"Wrote {PLOT_PREFIX}.png")
    print(f"Wrote {PLOT_DATA}")
    print(f"Wrote {SUMMARY_TABLE}")
    print(f"Wrote {COUNTRY_SUMMARY_TABLE}")


if __name__ == "__main__":
    main()
