#!/bin/bash

set -e

# Set a name for the output
for x in *.faa.gz; do n=\$x; done

# Make the MMSeqs2 database
mmseqs createdb *.faa.gz db

output_name="\$n.${min_identity}.${round}"

# Cluster the protein sequences
mmseqs cluster db \$output_name.cluster ./ \
    --min-seq-id ${min_identity / 100} \
    --max-seqs 100000 \
    -c ${min_coverage / 100}

# Get the representative sequences
mmseqs result2repseq db \$output_name.cluster \$output_name.rep
mmseqs result2flat db db \$output_name.rep \$output_name.dedup.fasta --use-fasta-header

gzip \$output_name.dedup.fasta
