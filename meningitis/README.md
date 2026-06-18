# Meningitis outbreak analysis

Methods for analysing the 2026 meningococcal outbreak in the United Kingdom.

All commands are included here, together with the most important files for reproducibility.

Versions used:
* atb command line tool 0.17.1 (https://github.com/AllTheBacteria/atb-cli)
* AllTheBacteria genomes releases r0.2, incremental 202408, plus incremental 202505
* Gubbins 3.4.3 (https://github.com/nickjcroucher/gubbins)
* LexicMap 0.9.0 (https://github.com/shenwei356/LexicMap)
* SeqKit 2.13.0 (https://github.com/shenwei356/seqkit)
* faqt 0.5.0 (https://github.com/martinghunt/faqt)


Download the reference genome:
```bash
faqt download-genome -o GCF_000191525.1.gff GCF_000191525.1
faqt to-fasta -i GCF_000191525.1.gff -o GCF_000191525.1.fa
```

Download all ST-485 isolates (script needs rauth installed):
```bash
python st_485_dl_auth.py \
    --key-name my_pubmlst_key \
    --output-dir genomes \
    --st 485 \
    --no-metadata
```

Get just the genomes we want:
```
cp genomes/neisseria_ST485_190673.fasta 190673_1926231.fas
cp genomes/neisseria_ST485_190682.fasta 190682_1930357.fas
cp genomes/neisseria_ST485_190683.fasta 190683_1930355.fas
cp genomes/neisseria_ST485_190684.fasta 190684_1930356.fas
cp genomes/neisseria_ST485_190685.fasta 190685_1930358.fas
```

Find the nearest genomes to the five outbreak genomes, using the `atb` command line wrapper around sketchlib:

```console
$ time atb sketch query --format tsv  --knn 500 1906*fas > atb_matches.tsv
Sketching 5 input file(s)...
🧬🖋️ sketchlib done in 0s
Querying ATB database (2,440,377 genomes)...
500 match(es) found.


real	0m7.406s
user	0m6.942s
sys	0m2.949s
```

We use `--knn 500` to get the top 500 hits, but the unique list of matching AllTheBacteria samples is shorter because the queries can match the same AllTheBacteria samples. (With `--knn 100` we actually get back 21 unique AllTheBacteria genomes.)

Get the sample accessions (deduplicate repeats):

```bash
awk 'NR>1 {print $2}' atb_matches.tsv | sort | uniq > atb_matches.accessions.txt
```

This is 100 matches:

```console
$ wc -l atb_matches.accessions.txt
     100 atb_matches.accessions.txt
```

Get the 100 matching genomes using `atb` commands:
```
atb query --sample-file atb_matches.accessions.txt > atb_matches.query_results.tsv
atb download --output-dir Matching_genomes --from atb_matches.query_results.tsv
```

Make a gubbins container definition file:

```console
cat > gubbins_3.4.3.def <<'EOF'
Bootstrap: docker
From: condaforge/miniforge3:latest

%post
    mamba create -y -n gubbins -c conda-forge -c bioconda python=3.10 gubbins=3.4.3
    mamba clean -afy

%environment
    export PATH=/opt/conda/envs/gubbins/bin:/opt/conda/bin:$PATH

%runscript
    exec /opt/conda/envs/gubbins/bin/run_gubbins.py "$@"
EOF
```

Build the container:
```bash
sudo singularity build gubbins_3.4.3.sif gubbins_3.4.3.def
```


Make gubbins input list file:
```bash
ls 1906*.fas | awk -F. '{print $1"\t"$0}' > gubbins.input.list
awk '{print $1"\tMatching_genomes/"$1".fa.gz"}' atb_matches.accessions.txt >> gubbins.input.list
```

Run gubbins:
```bash
singularity exec gubbins_3.4.3.sif generate_ska_alignment.py --reference GCF_000191525.1.fa --input gubbins.input.list --out gubbins.out.aln
singularity exec gubbins_3.4.3.sif run_gubbins.py --filter-percentage 30 --prefix gubbins.out gubbins.out.aln 2>&1 | tee gubbins.run.log
```

View these files in Phandango: `gubbins.out.recombination_predictions.gff`, `GCF_000191525.1.gff`, `gubbins.out.final_tree.tre`.

Get the nucleotide sequences of the two genes of interest:

```bash
seqkit subseq -r 2131760:2132752 GCF_000191525.1.fa > WP_002248791.1.fa
seqkit subseq -r 638007:638396 GCF_000191525.1.fa > WP_002226625.1.fa
```

Run lexicmap (using lexicmap database built from releases 0.2 plus 2024-08 plus 2025-05):

```bash
lexicmap search -d atb-202505.lmi -j 4 -o WP_002226625.1.lexicmap.tsv.gz WP_002226625.1.fa
lexicmap search -d atb-202505.lmi -j 4 -o WP_002248791.1.lexicmap.tsv.gz WP_002248791.1.fa
```


Gather basic stats and write file of sample accessions for querying in AllTheBacteria to get species and dates:

```bash
for x in WP_002226625.1.lexicmap.tsv.gz WP_002248791.1.lexicmap.tsv.gz ; do
    echo "---------  $x  ---------";
    a=`gunzip -c $x | wc -l`
    echo "total matches: $((a-1))"
    a=`gunzip -c $x | awk '{print $4}' | sort -u | wc -l`
    echo "total matching samples: $((a-1))"
    gunzip -c $x | awk -F"\t" '$6==100 && $11==100 {print $4}' > $x.perfect_match_samples.txt
    a=`wc -l $x.perfect_match_samples.txt | awk '{print $1}'`
    echo "Samples with perfect match: $a"
done
```

Output:
```console
---------  WP_002226625.1.lexicmap.tsv.gz  ---------
total matches: 104011
total matching samples: 103694
Samples with perfect match: 2020
---------  WP_002248791.1.lexicmap.tsv.gz  ---------
total matches: 169481
total matching samples: 101646
Samples with perfect match: 4636
```

Query AllTheBacteria:

```bash
atb query --columns sample_accession,hq_filter,sylph_species,country,collection_date --sample-file WP_002226625.1.lexicmap.tsv.gz.perfect_match_samples.txt > WP_002226625.1.lexicmap.tsv.gz.perfect_match_samples.atb_meta.tsv

atb query --columns sample_accession,hq_filter,sylph_species,country,collection_date --sample-file WP_002248791.1.lexicmap.tsv.gz.perfect_match_samples.txt > WP_002248791.1.lexicmap.tsv.gz.perfect_match_samples.atb_meta.tsv
```

Make the plot:
```console
$ python plot_allele_country_dates.py
Per-allele plotted counts:
WP_002226625.1: 239 genomes, 13 countries, 1967-2023
WP_002248791.1: 1192 genomes, 16 countries, 1978-2025
Wrote supplementary_allele_country_dates.pdf
Wrote supplementary_allele_country_dates.png
Wrote supplementary_allele_country_year_counts.tsv
Wrote supplementary_allele_country_date_summary.tsv
Wrote supplementary_allele_country_summary.tsv
```

