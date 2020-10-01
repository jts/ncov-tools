#
# Wrapper for accessing the config, with sensible defaults
#

#
def get_bam_pattern():
    return config.get("bam_pattern", "{data_root}/{sample}.sorted.bam")

#
def get_consensus_pattern():
    if config['platform'] == 'illumina':
        return config.get("consensus_pattern", "{data_root}/{sample}.primertrimmed.consensus.fa")
    elif config['platform'] == 'oxford-nanopore':
        return config.get("consensus_pattern", "{data_root}/{sample}.consensus.fasta")

#
def get_variants_pattern():
    if config['platform'] == 'illumina':
        return config.get("variants_pattern", "{data_root}/{sample}.variants.tsv")
    elif config['platform'] == 'oxford-nanopore':
        return config.get("variants_pattern", "{data_root}/{sample}.pass.vcf.gz")

#
def get_metadata_file(wildcards):
    return config.get("metadata", "").format(data_root=config["data_root"])

#
def get_run_name(wildcards=None):
    return config.get("run_name", "default")

# Get the path to the BED file containing amplicon primers.
def get_primer_bed(wildcards):
    return config.get("primer_bed", "")
    
# Get the number of bases to offset the primer when removing primers from amplicon BED
def get_primer_offset(wildcards):
    return config.get("primer_offset", "0")

# Get option bed_type from config.yaml
def get_primer_bed_type_opt(wildcards):
    return config.get("bed_type", "unique_amplicons")

#
def get_completeness_threshold(wildcards):
    return config.get("completeness_threshold", 0.75)

#
def get_sarscov2db_opt(wildcards):
    return config.get("sarscov2db", "")
