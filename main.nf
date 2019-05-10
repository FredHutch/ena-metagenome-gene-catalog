#!/usr/bin/env nextflow

params.biome_name = "Host-associated"
params.min_coverage = 50
params.output_folder = "./"
min_identity_ch = Channel.from( 100, 90, 80, 70, 60, 50 )

// Fetch the list of host-associated studies from ENA
process fetchHostAssociatedStudies {
  container "quay.io/fhcrc-microbiome/python-pandas@sha256:39993ba37c44368d1a5752cf6b96f8172e69eb109374722bd6914c29a79565c6"
  cpus 2
  memory "4 GB"
  
  input:
  val biome_name from params.biome_name
  
  output:
  file "list_of_assemblies.txt" into list_of_assemblies

  """
#!/usr/bin/python3

import requests


def get(url):
    r = requests.get(url)
    return r.json()


def parse_json(d, key_list):
    for k in key_list:
        assert k in d
        d = d[k]
    return d


def get_all(url, item=["data"], next_key=["links", "next"], n=None):
    ix = 1
    total_set = []
    d = get(url)
    total_set.extend(parse_json(d, item))

    while parse_json(d, next_key) is not None:
        if n is not None and ix == n:
            break
        d = get(parse_json(d, next_key))
        total_set.extend(parse_json(d, item))
        ix += 1

    return total_set


# GET ALL HOST ASSOCIATED STUDIES
assembly_list = get_all(
    "https://www.ebi.ac.uk/metagenomics/api/v1/assemblies?lineage=root%3A${biome_name}&page=1"
)

print("Fetched a total of " + str(len(assembly_list)) + " ${biome_name} assemblies")

open("list_of_assemblies.txt", "wt").write("\\n".join([
    parse_json(assembly, ["relationships", "analyses", "links", "related"])
    for assembly in assembly_list
]))

  """

}


// Break up the list of assemblies
process splitAssemblyList {
    container "ubuntu:16.04"
    cpus 1
    memory "1 GB"

    input:
    file list_of_assemblies

    output:
    file "*" into groups_of_assemblies

    """
    split -l 10 ${list_of_assemblies}
    """
}


// Download the CDS per assembly
process fetchCDS {
  container "quay.io/fhcrc-microbiome/python-pandas@sha256:39993ba37c44368d1a5752cf6b96f8172e69eb109374722bd6914c29a79565c6"
  cpus 1
  memory "1 GB"
  
  input:
  file assembly_url_list from groups_of_assemblies.flatten()
  
  output:
  file "*.faa.gz" into cds_ch

  """
#!/usr/bin/python3

import os
import requests


def get(url):
    r = requests.get(url)
    return r.json()


def parse_json(d, key_list):
    for k in key_list:
        assert k in d, d.keys()
        d = d[k]
    return d


def get_all(url, item=["data"], next_key=["links", "next"], n=None):
    ix = 1
    total_set = []
    d = get(url)
    total_set.extend(parse_json(d, item))

    while parse_json(d, next_key) is not None:
        if n is not None and ix == n:
            break
        d = get(parse_json(d, next_key))
        total_set.extend(parse_json(d, item))
        ix += 1

    return total_set

def download_file(file_url):
    print("Downloading " + file_url)
    filename = file_url.split("/")[-1]
    assert os.path.exists(filename) is False
    r = requests.get(file_url, allow_redirects=True)
    with open(filename, 'wb') as fo:
        fo.write(r.content)

for assembly_url in open("${assembly_url_list}").readlines():
    assembly_url = assembly_url.strip()
    print("Processing " + assembly_url)

    # Get available analyses
    break_out = False
    for analysis in get_all(assembly_url):
        if break_out:
            break
        download_url = parse_json(analysis, ["relationships", "downloads", "links", "related"])

        for download in get_all(download_url):
            if parse_json(download, ["attributes", "description", "label"]) == "Predicted CDS without annotation":
                faa_url = parse_json(download, ["links", "self"])
                download_file(faa_url)
                if break_out:
                    break

  """

}


// Combine 100% identical sequences
process deduplicateCDS {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    cpus 2
    memory "4 GB"

    input:
    file "*" from cds_ch
    val min_identity from 100
    val min_coverage from 50
    
    output:
    file "*.dedup.fasta.gz" into dedup_ch

    """
#!/bin/bash

set -e

cat *.faa.gz > batch_cds.faa.gz
for x in *.faa.gz; do n=\$x; done

# Make the MMSeqs2 database
mmseqs createdb batch_cds.faa.gz db

# Cluster the protein sequences
mmseqs cluster db \$n.${min_identity}.cluster ./ \
    --min-seq-id ${min_identity / 100} \
    --max-seqs 100000 \
    -c ${min_coverage / 100}

# Get the representative sequences
mmseqs result2repseq db \$n.${min_identity}.cluster \$n.${min_identity}.rep
mmseqs result2flat db db \$n.${min_identity}.rep \$n.${min_identity}.dedup.fasta --use-fasta-header

gzip \$n.${min_identity}.dedup.fasta

    """
}


process combineCDS {
    container "ubuntu:16.04"
    cpus 1
    memory "8 GB"
    
    input:
    file "*" from dedup_ch.collect()
    
    output:
    file "all_CDS.faa.gz" into all_cds

    """
set -e

cat *faa.gz > all_CDS.faa.gz
gzip -t all_CDS.faa.gz
    """
}

process clusterCDS {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    cpus 4
    memory "8 GB"
    publishDir "${params.output_folder}"
    
    input:
    file all_cds
    val min_identity from min_identity_ch
    val min_coverage from params.min_coverage
    
    output:
    file "mmseqs.${min_identity}.tsv.gz"
    file "mmseqs.${min_identity}.rep.fasta.gz"

    """
#!/bin/bash

set -e

# Make the MMSeqs2 database
mmseqs createdb ${all_cds} db

# Cluster the protein sequences
mmseqs cluster db mmseqs.${min_identity}.cluster ./ \
    --min-seq-id ${min_identity / 100} \
    --max-seqs 100000 \
    -c ${min_coverage / 100}

# Make TSV output for clustering
mmseqs createtsv db db mmseqs.${min_identity}.cluster mmseqs.${min_identity}.tsv

# Get the representative sequences
mmseqs result2repseq db mmseqs.${min_identity}.cluster mmseqs.${min_identity}.rep
mmseqs result2flat db db mmseqs.${min_identity}.rep mmseqs.${min_identity}.rep.fasta --use-fasta-header

gzip mmseqs.${min_identity}.tsv
gzip mmseqs.${min_identity}.rep.fasta

    """
}