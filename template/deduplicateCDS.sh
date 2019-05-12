#!/bin/bash

set -e

# Set a name for the output
for x in *.faa.gz; do n=\$x; done

# Make the MMSeqs2 database
mmseqs createdb *.faa.gz db

# Cluster the protein sequences
mmseqs cluster db \$n.${min_identity}.cluster ./ \
    --min-seq-id ${min_identity / 100} \
    --max-seqs 100000 \
    -c ${min_coverage / 100}

# Get the representative sequences
mmseqs result2repseq db \$n.${min_identity}.cluster \$n.${min_identity}.rep
mmseqs result2flat db db \$n.${min_identity}.rep \$n.${min_identity}.dedup.fasta --use-fasta-header

gzip \$n.${min_identity}.dedup.fasta
