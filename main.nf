#!/usr/bin/env nextflow

params.biome_name = "Host-associated"
params.min_coverage = 50
params.output_folder = "./"
min_identity_ch = Channel.from( 99, 90, 80, 70, 60, 50 )

// Fetch the number of pages of host-associated studies from ENA
process pagesHostAssociatedStudies {
  container "quay.io/fhcrc-microbiome/python-pandas@sha256:39993ba37c44368d1a5752cf6b96f8172e69eb109374722bd6914c29a79565c6"
  cpus 2
  memory "4 GB"
  
  input:
  val biome_name from params.biome_name
  
  output:
  file "page.*.txt" into list_of_urls

  """
#!/usr/bin/python3

import requests
import time


def get(url, attempts=5, delay=60):
    i = 0
    while i < attempts:
        try:
            r = requests.get(url)
            return r.json()
        except:
            time.sleep(delay)
        i += 1


def parse_json(d, key_list):
    for k in key_list:
        assert k in d
        d = d[k]
    return d


# Get the first page of assemblies
d = get(
    "https://www.ebi.ac.uk/metagenomics/api/v1/assemblies?lineage=root%3A${biome_name}"
)

# Figure out how many pages there are
n_pages = parse_json(d, ["meta", "pagination", "pages"])

print("There are a total of " + str(n_pages) + " pages of ${biome_name} assemblies")

for ix in range(1, n_pages + 1):
    with open("page." + str(ix) + ".txt", "wt") as fo:
        fo.write("https://www.ebi.ac.uk/metagenomics/api/v1/assemblies?lineage=root%3A${biome_name}&page=" + str(ix))

  """

}


// Fetch the list of host-associated studies from ENA
process fetchHostAssociatedStudies {
  container "quay.io/fhcrc-microbiome/python-pandas@sha256:39993ba37c44368d1a5752cf6b96f8172e69eb109374722bd6914c29a79565c6"
  cpus 1
  memory "1 GB"
  errorStrategy 'retry'
  
  input:
  file url from list_of_urls.flatten()
  val biome_name from params.biome_name

  output:
  file "assemblies*" into groups_of_assemblies

  """
#!/usr/bin/python3

import requests
import time


def get(url, attempts=5, delay=60):
    i = 0
    while i < attempts:
        try:
            r = requests.get(url)
            return r.json()
        except:
            time.sleep(delay)
        i += 1


def parse_json(d, key_list):
    for k in key_list:
        assert k in d, (k, d.keys())
        d = d[k]
    return d


# Get the page to read in this process
with open("${url}", "rt") as f:
    page_url = f.readline().strip()

# Get the host-associated assemblies from this page
assembly_list = get(page_url)

print("Fetched " + str(len(assembly_list)) + " ${biome_name} assemblies")

open("assemblies.${url}", "wt").write("\\n".join([
    parse_json(assembly, ["relationships", "analyses", "links", "related"])
    for assembly in assembly_list["data"]
]))

  """

}


// Download the CDS per assembly
process fetchCDS {
  container "quay.io/fhcrc-microbiome/python-pandas@sha256:39993ba37c44368d1a5752cf6b96f8172e69eb109374722bd6914c29a79565c6"
  cpus 1
  memory "8 GB"
  errorStrategy 'retry'
  
  input:
  file assembly_url_list from groups_of_assemblies.flatten()
  
  output:
  file "*.faa.gz" into dedup_one

  afterScript "rm -r *"

  """
#!/usr/bin/python3

import os
import requests
import time


def get(url, attempts=5, delay=60):
    i = 0
    while i < attempts:
        try:
            r = requests.get(url)
            return r.json()
        except:
            time.sleep(delay)
        i += 1


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

def download_file(file_url, filename):
    print("Downloading " + file_url)
    print("Local destination " + filename)
    assert os.path.exists(filename) is False
    r = requests.get(file_url, allow_redirects=True)
    with open(filename, 'wb') as fo:
        fo.write(r.content)

all_filenames = set([])

for assembly_url in open("${assembly_url_list}").readlines():
    assembly_url = assembly_url.strip()
    print("Processing " + assembly_url)

    # Get available analyses
    for analysis in get_all(assembly_url):

        analysis_id = parse_json(analysis, ["id"])

        download_url = parse_json(analysis, ["relationships", "downloads", "links", "related"])

        for download in get_all(download_url):
            if parse_json(download, ["attributes", "description", "label"]) == "Predicted CDS without annotation":
                faa_url = parse_json(download, ["links", "self"])

                # Format the filename by the analysis' unique ID
                filename = analysis_id + ".CDS.unannotated.faa.gz"

                # Make sure this file hasn't been downloaded yet
                if filename not in all_filenames:

                    download_file(faa_url, filename)

                    all_filenames.add(filename)

  """

}


// Combine 100% identical sequences in multiple rounds
process deduplicateRoundOne {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    cpus 16
    memory "120 GB"
    errorStrategy 'retry'

    input:
    file "*" from dedup_one.flatten().collate(5)
    val min_identity from 99
    val min_coverage from 50
    val round from 1
    
    output:
    file "*.dedup.fasta.gz" into dedup_two

    afterScript "rm -r *"

    script:
    template "deduplicateCDS.sh"
}


process deduplicateRoundTwo {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    cpus 16
    memory "120 GB"
    errorStrategy 'retry'

    input:
    file "*" from dedup_two.flatten().collate(5)
    val min_identity from 99
    val min_coverage from 50
    val round from 2
    
    output:
    file "*.dedup.fasta.gz" into dedup_three

    afterScript "rm -r *"

    script:
    template "deduplicateCDS.sh"
}


process deduplicateRoundThree {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    cpus 16
    memory "120 GB"
    errorStrategy 'retry'

    input:
    file "*" from dedup_three.flatten().collate(2)
    val min_identity from 99
    val min_coverage from 50
    val round from 3
    
    output:
    file "*.dedup.fasta.gz" into dedup_ch

    afterScript "rm -r *"

    script:
    template "deduplicateCDS.sh"
}


process combineCDS {
    container "ubuntu:16.04"
    cpus 16
    memory "120 GB"
    errorStrategy 'retry'
    
    input:
    file "*" from dedup_ch.collect()
    
    output:
    file "all_CDS.faa.gz" into all_cds

    afterScript "rm -r *"

    """
set -e

cat *faa.gz > all_CDS.faa.gz
gzip -t all_CDS.faa.gz
    """
}

process clusterCDS {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    cpus 16
    memory "120 GB"
    publishDir "${params.output_folder}"
    errorStrategy 'retry'
    
    input:
    file all_cds
    val min_identity from min_identity_ch
    val min_coverage from params.min_coverage
    
    output:
    file "mmseqs.${min_identity}.tsv.gz"
    file "mmseqs.${min_identity}.rep.fasta.gz"

    afterScript "rm -r *"

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