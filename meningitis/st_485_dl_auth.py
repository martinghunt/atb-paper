#!/usr/bin/env python3
"""
Download all Neisseria ST-485 genomes from pubMLST as FASTA files, with OAuth1
authentication via BIGSdb/PubMLST credentials.

Usage:
    python st_485_dl_auth.py --key-name <name> [--output-dir ./genomes] [--delay 0.5]
"""

from __future__ import annotations

import argparse
import configparser
import json
import logging
import os
import stat
import time
from pathlib import Path

from rauth import OAuth1Service, OAuth1Session

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

BASE_WEB = {
    "PubMLST": "https://pubmlst.org/bigsdb",
    "Pasteur": "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl",
}
BASE_API = {
    "PubMLST": "https://rest.pubmlst.org",
    "Pasteur": "https://bigsdb.pasteur.fr/api",
}

DEFAULT_DB = "pubmlst_neisseria_isolates"
DEFAULT_SCHEME_ID = 1
DEFAULT_ST = 485
DEFAULT_PAGE_SIZE = 100


def check_dir(directory: Path) -> None:
    if directory.is_dir():
        if os.access(directory, os.W_OK):
            return
        raise PermissionError(f"The token directory '{directory}' exists but is not writable.")
    directory.mkdir(parents=True, exist_ok=True)
    os.chmod(directory, stat.S_IRWXU)


def read_token(token_dir: Path, token_type: str, key_name: str) -> tuple[str | None, str | None]:
    file_path = token_dir / f"{token_type}_tokens"
    if file_path.is_file():
        config = configparser.ConfigParser(interpolation=None)
        config.read(file_path)
        if config.has_section(key_name):
            return config[key_name]["token"], config[key_name]["secret"]
    return None, None


def write_token(token_dir: Path, token_type: str, key_name: str, token: str, secret: str) -> None:
    file_path = token_dir / f"{token_type}_tokens"
    config = configparser.ConfigParser(interpolation=None)
    if file_path.is_file():
        config.read(file_path)
    config[key_name] = {"token": token, "secret": secret}
    with open(file_path, "w") as configfile:
        config.write(configfile)


def get_client_credentials(token_dir: Path, key_name: str) -> tuple[str, str]:
    config = configparser.ConfigParser(interpolation=None)
    file_path = token_dir / "client_credentials"
    client_id = None
    if file_path.is_file():
        config.read(file_path)
        if config.has_section(key_name):
            client_id = config[key_name]["client_id"]
            client_secret = config[key_name]["client_secret"]

    if not client_id:
        client_id = input("Enter client id: ").strip()
        while len(client_id) != 24:
            print("Client ids are exactly 24 characters long.")
            client_id = input("Enter client id: ").strip()
        client_secret = input("Enter client secret: ").strip()
        while len(client_secret) != 42:
            print("Client secrets are exactly 42 characters long.")
            client_secret = input("Enter client secret: ").strip()
        config[key_name] = {"client_id": client_id, "client_secret": client_secret}
        with open(file_path, "w") as configfile:
            config.write(configfile)

    return client_id, client_secret


def get_service(site: str, db: str, client_key: str, client_secret: str) -> OAuth1Service:
    request_token_url = f"{BASE_API[site]}/db/{db}/oauth/get_request_token"
    access_token_url = f"{BASE_API[site]}/db/{db}/oauth/get_access_token"
    return OAuth1Service(
        name="bigsdb_downloader_wrapper",
        consumer_key=client_key,
        consumer_secret=client_secret,
        request_token_url=request_token_url,
        access_token_url=access_token_url,
        base_url=BASE_API[site],
    )


def get_new_access_token(
    site: str, db: str, token_dir: Path, key_name: str, client_key: str, client_secret: str
) -> tuple[str, str]:
    service = get_service(site, db, client_key, client_secret)
    r = service.get_raw_request_token(
        params={"oauth_callback": "oob"}, headers={"User-Agent": "BIGSdb downloader wrapper"}
    )
    if r.status_code != 200:
        raise RuntimeError("Failed to get new request token.")

    request_token = r.json()["oauth_token"]
    request_secret = r.json()["oauth_token_secret"]

    print(
        "Please log in using your user account at "
        f"{BASE_WEB[site]}?db={db}&page=authorizeClient&oauth_token={request_token} "
        "using a web browser to obtain a verification code."
    )
    verifier = input("Please enter verification code: ").strip()
    r = service.get_raw_access_token(
        request_token,
        request_secret,
        params={"oauth_verifier": verifier},
        headers={"User-Agent": "BIGSdb downloader wrapper"},
    )
    if r.status_code != 200:
        raise RuntimeError("Failed to get new access token.")

    token = r.json()["oauth_token"]
    secret = r.json()["oauth_token_secret"]
    write_token(token_dir, "access", key_name, token, secret)
    return token, secret


def get_session_token(
    site: str, db: str, token_dir: Path, key_name: str
) -> tuple[str, str]:
    client_key, client_secret = get_client_credentials(token_dir, key_name)
    access_token, access_secret = read_token(token_dir, "access", key_name)
    if not access_token or not access_secret:
        access_token, access_secret = get_new_access_token(
            site, db, token_dir, key_name, client_key, client_secret
        )

    url = f"{BASE_API[site]}/db/{db}/oauth/get_session_token"
    session_request = OAuth1Session(
        client_key,
        client_secret,
        access_token=access_token,
        access_token_secret=access_secret,
    )
    r = session_request.get(url, headers={"User-Agent": "BIGSdb downloader wrapper"})
    if r.status_code != 200:
        raise RuntimeError(f"Failed to get new session token: {r.text}")

    token = r.json()["oauth_token"]
    secret = r.json()["oauth_token_secret"]
    write_token(token_dir, "session", key_name, token, secret)
    return token, secret


def get_oauth_session(site: str, db: str, token_dir: Path, key_name: str) -> OAuth1Session:
    client_key, client_secret = get_client_credentials(token_dir, key_name)
    token, secret = read_token(token_dir, "session", key_name)
    if not token or not secret:
        token, secret = get_session_token(site, db, token_dir, key_name)
    return OAuth1Session(
        client_key,
        client_secret,
        access_token=token,
        access_token_secret=secret,
    )


def get_isolate_uris(
    session: OAuth1Session, site: str, db: str, scheme_id: int, st: int, page_size: int
) -> list[str]:
    url = f"{BASE_API[site]}/db/{db}/isolates/search"
    payload = {f"scheme.{scheme_id}.ST": st}
    uris: list[str] = []
    page = 1

    while True:
        r = session.post(
            url,
            params={"page_size": page_size, "page": page},
            data=json.dumps(payload),
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
                "User-Agent": "BIGSdb downloader wrapper",
            },
            header_auth=True,
        )
        if r.status_code != 200 and r.status_code != 201:
            raise RuntimeError(f"Search failed: {r.text}")
        data = r.json()

        if page == 1:
            log.info("Query matched %s isolates", data.get("records", "?"))

        batch = data.get("isolates", [])
        uris.extend(batch)
        log.info("Page %d: %d URIs collected so far", page, len(uris))

        if len(batch) < page_size:
            break
        page += 1

    return uris


def download_fasta(
    session: OAuth1Session, site: str, db: str, isolate_id: str, dest: Path
) -> None:
    url = f"{BASE_API[site]}/db/{db}/isolates/{isolate_id}/contigs_fasta"
    r = session.get(
        url,
        params={"header": "original_designation"},
        headers={"User-Agent": "BIGSdb downloader wrapper"},
    )
    if r.status_code != 200:
        raise RuntimeError(f"Download failed for isolate {isolate_id}: {r.text}")
    dest.write_text(r.text)


def download_metadata(
    session: OAuth1Session, site: str, db: str, isolate_id: str, dest: Path
) -> None:
    url = f"{BASE_API[site]}/db/{db}/isolates/{isolate_id}"
    r = session.get(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "BIGSdb downloader wrapper",
        },
    )
    if r.status_code != 200:
        raise RuntimeError(f"Metadata fetch failed for isolate {isolate_id}: {r.text}")
    dest.write_text(json.dumps(r.json()))


def main(
    site: str,
    db: str,
    scheme_id: int,
    st: int,
    output_dir: Path,
    delay: float,
    page_size: int,
    key_name: str,
    token_dir: Path,
    metadata: bool,
    fasta: bool,
) -> None:
    check_dir(token_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    session = get_oauth_session(site, db, token_dir, key_name)

    log.info("Querying ST-%d isolates from %s", st, db)
    isolate_uris = get_isolate_uris(session, site, db, scheme_id, st, page_size)
    log.info("Found %d isolates with ST-%d", len(isolate_uris), st)

    downloaded, skipped = 0, 0
    for uri in isolate_uris:
        isolate_id = uri.rstrip("/").split("/")[-1]
        dest = output_dir / f"neisseria_ST{st}_{isolate_id}.fasta"
        meta_dest = output_dir / f"neisseria_ST{st}_{isolate_id}.json"

        if not fasta and not metadata:
            log.debug("Skipping %s (nothing requested)", isolate_id)
            skipped += 1
            continue

        try:
            did_work = False
            if fasta:
                if dest.exists():
                    log.debug("Skipping %s (already exists)", dest.name)
                else:
                    download_fasta(session, site, db, isolate_id, dest)
                    did_work = True
            if metadata:
                if meta_dest.exists():
                    log.debug("Skipping %s (already exists)", meta_dest.name)
                else:
                    download_metadata(session, site, db, isolate_id, meta_dest)
                    did_work = True
            log.info("Downloaded isolate %s -> %s", isolate_id, dest.name)
            if did_work:
                downloaded += 1
            else:
                skipped += 1
        except Exception as e:
            log.error("Error for isolate %s: %s", isolate_id, e)
            skipped += 1

        if delay > 0:
            time.sleep(delay)

    log.info("Done. Downloaded: %d, Skipped/failed: %d", downloaded, skipped)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download MLST isolate FASTA files from PubMLST with OAuth1 auth"
    )
    parser.add_argument("--site", choices=["PubMLST", "Pasteur"], default="PubMLST")
    parser.add_argument("--db", default=DEFAULT_DB)
    parser.add_argument("--scheme-id", type=int, default=DEFAULT_SCHEME_ID)
    parser.add_argument("--st", type=int, default=DEFAULT_ST)
    parser.add_argument("--page-size", type=int, default=DEFAULT_PAGE_SIZE)
    parser.add_argument("--output-dir", type=Path, default=Path("./genomes"))
    parser.add_argument("--delay", type=float, default=0.5)
    parser.add_argument("--key-name", required=True, help="API key name used in .bigsdb_tokens")
    parser.add_argument("--token-dir", type=Path, default=Path("./.bigsdb_tokens"))
    meta_group = parser.add_mutually_exclusive_group()
    meta_group.add_argument("--metadata", dest="metadata", action="store_true", help="Also save isolate metadata JSON")
    meta_group.add_argument("--no-metadata", dest="metadata", action="store_false", help="Skip metadata JSON download")
    parser.set_defaults(metadata=True)
    fasta_group = parser.add_mutually_exclusive_group()
    fasta_group.add_argument("--fasta", dest="fasta", action="store_true", help="Download FASTA sequences")
    fasta_group.add_argument("--no-fasta", dest="fasta", action="store_false", help="Skip FASTA download")
    parser.set_defaults(fasta=True)
    args = parser.parse_args()

    main(
        args.site,
        args.db,
        args.scheme_id,
        args.st,
        args.output_dir,
        args.delay,
        args.page_size,
        args.key_name,
        args.token_dir,
        args.metadata,
        args.fasta,
    )
