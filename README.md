# ncov-tools

Various tools and QC plots for working with coronavirus sequencing data.

## Installation

Download the package:
```
git clone git@github.com:jts/ncov-tools.git
cd ncov-tools
```

To use this package, install the dependencies using conda:
```
conda env create -f environment.yml
```

Alternatively, if install times are very slow using conda, we recommend using
the conda wrapper: [mamba](https://github.com/TheSnakePit/mamba).

Install mamba as follows:
```
conda install -c conda-forge mamba
```

Then create the ncov-tools environment using mamba

```
mamba env create -f environment.yml
```

Either way, if you used conda directly or mamba, activate 
the conda package and install the `parser` package:
```
conda activate ncov-qc
cd parser
pip install -r requirements.txt
pip install .
```

## Configuration

This package is implemented as a snakemake pipeline, so requires a `config.yaml` file to describe where the input files are. To generate QC plots, a bam file with reads mapped to a reference  genome is required. Consensus sequences (FASTA) are needed to generate a phylogenetic tree with associated mutations.

As an example, let's say your data is laid out in the following structure:

```
   run_200430/
     sampleA.sorted.bam
     sampleA.consensus.fasta
     sampleB.sorted.bam
     sampleB.consensus.fasta
   resources/
     artic_reference.fasta
     artic_amplicons.bed
```

Then your config.yaml should look like:

```
# path to the top-level directory containing the analysis results
data_root: run_200430

# optionally the plots can have a "run name" prefix. If this is not defined the prefix will be "default"
run_name: my_run

# path to the file containing the amplicon regions (not the primer sites, the actual amplicons)
amplicon_bed: resources/artic_amplicons.bed

# path to the nCov reference genome
reference_genome: resources/artic_reference.fasta

# the naming convention for the bam files
# this can use the variables {data_root} (as above) and {sample}
# As per the example above, this will expand to run_200430/sampleA.sorted.bam for sampleA
bam_pattern: "{data_root}/{sample}.sorted.bam"

# the naming convention for the consensus sequences
consensus_pattern: "{data_root}/{sample}.consensus.fasta"

# when building a tree of the consensus genomes you can optionally include other sequences
# in the tree by providing a fasta file here
tree_include_consensus: some_genomes_from_gisaid.fasta

# some plots can by annotated with external metadata. this
# file contains the metadata in simple tab-delimited format
# one column must be 'sample'
metadata: metadata.tsv

#
# set this flag to true to include lineage assignments with pangolin in the output plots
#
assign_lineages: true

#
# this flag should identify the sequencing platform used, valid options are:
# `illumina` or `oxford-nanopore`
#
platform: illumina
```

## Running

After configuration, you can run the pipeline using Snakemake

```
# Build the sequencing QC plots (coverage, allele frequencies)
snakemake -s qc/Snakefile all_qc_sequencing

# Build the analysis QC plots (tree with annotated mutations)
snakemake -s qc/Snakefile all_qc_analysis
```

The results will be placed in the `plots` directory.

To generate a `tsv` file of QC metrics per sample:
```
snakemake -s qc/Snakefile all_qc_summary --cores 1
```

## Credit and Acknowledgements

The tree-with-SNPs plot was inspired by a plot shared by Mads Albertsen.

