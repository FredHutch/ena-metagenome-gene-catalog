# ENA Metagenome Gene Catalog

Download metagenomic assemblies from ENA and assemble a gene catalog

### Workflow

  1. Identify all of the host-associated metagenomes with the ENA API
  2. Download each metagenome assembly and extract the CDS (gene sequences)
  3. Aggregate all gene sequences
  4. Cluster genes at a few different levels of amino acid identity
  5. Return the non-redundant gene catalogs

### Output

For each level of amino acid identity (e.g. 100%, 90%, 80%):

  * FASTA of protein sequences
  * TSV linking each input gene to the final cluster
