#
# Wrapper for accessing the config, with sensible defaults
#

import re

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
    else:
        sys.stderr.write("Error, unknown platform %s\n" % (config['platform']))
        sys.exit(1)
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

# get the primer name prefix from the config.yaml file
def get_primer_prefix(wildcards):
    return config.get("primer_prefix", "nCoV-2019")

def get_snp_tree_flag(wildcards=None):
    return config.get("build_snp_tree", True)

#
def get_completeness_threshold(wildcards):
    return config.get("completeness_threshold", 0.75)

#
def get_sarscov2db_opt(wildcards):
    return config.get("sarscov2db", "")

#
def get_watch_mutation_set(wildcards):
    return config.get("mutation_set", "spike_mutations")

#
def get_pangolin_analysis_mode(wildcards):
    return config.get("pango_analysis_mode", "accurate")

def get_pangolin_analysis_mode_opt(wildcards):
    analysis_mode = get_pangolin_analysis_mode(wildcards)
    pangolin_version = config.get('pangolin_version', '')
    if pangolin_version == '4' or pangolin_version == '':
        return f'--analysis-mode {analysis_mode}'
    else:
        return ''

#
def get_voc_pango_lineages(wildcards):
    return config.get("voc_pango_lineages", "B.1.1.7,B.1.351,P.1,B.1.617.2")

#
def get_data_directory(wildcards=None):
    return config["data_root"]

#
def get_variant_format_pattern(wildcards=None):
    variant_format_pattern = config['variants_pattern']
    _pattern = re.sub('{data_root}/{sample}', '', variant_format_pattern)
    return _pattern

