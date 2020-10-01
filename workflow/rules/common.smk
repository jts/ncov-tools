#
# Helper functions and rules
#
include: "defaults.smk"

configfile: "config.yaml"

def get_sample_names():

    # if defined in the config, use that
    # otherwise we try to auto-detect based on the bam names
    if "samples" in config:
        return config["samples"]

    # get all bams in the data directory based on the pattern
    # the bam files follow
    pattern = get_bam_pattern()

    # form a glob we can use to get the bam files for each sample
    gs = pattern.format(data_root=config["data_root"], sample="*")

    bams = glob.glob(gs)
    samples = list()
    for b in bams:
        f = os.path.basename(b)
        fields = f.split(".")
        samples.append(fields[0])

    # remove any possible duplicates
    return list(set(samples))

def get_negative_control_samples():
    if "negative_control_samples" in config:
        return config["negative_control_samples"]
    else:
        return []

# negative controls may not have bams so we need to filter
# out those that don't have alignments
def get_negative_control_bed(wildcards):
    bp = get_bam_pattern()
    out = list()
    for s in get_negative_control_samples():
        bam = bp.format(data_root=config["data_root"], sample=s)
        if os.path.exists(bam):
            out.append("qc_sequencing/{sample}.amplicon_base_coverage.bed".format(sample=s))
    return out

def get_amplicon_bed(wildcards):
    return config["amplicon_bed"].format(data_root=config["data_root"])

def get_reference_genome(wildcards):
    return config["reference_genome"].format(data_root=config["data_root"])

def get_reference_genome_fai(wildcards):
    return get_reference_genome(wildcards) + ".fai"

def get_bam_for_sample(wildcards):
    bam_pattern = get_bam_pattern()
    return bam_pattern.format(data_root=config["data_root"], sample=wildcards.sample)

def get_primer_trimmed_bam_for_sample(wildcards):

    # If the pattern has been provided in the config, use it
    if "primer_trimmed_bam_pattern" in config:
        return config["primer_trimmed_bam_pattern"].format(data_root=config["data_root"], sample=wildcards.sample)

    # Otherwise, rely on what we expect the platform to produce
    if "platform" not in config:
        sys.stderr.write("Error: platform not defined in config")
        sys.exit(1)

    if config['platform'] == 'oxford-nanopore':
        return "%s/{sample}.primertrimmed.rg.sorted.bam" % config['data_root']
    elif config['platform'] == 'illumina':
        return "%s/{sample}.mapped.primertrimmed.sorted.bam" % config['data_root']
    else:
        sys.stderr.write("Error: unrecognized platform")
        sys.exit(1)

def get_completeness_threshold_opt(wildcards):
    completeness = get_completeness_threshold(wildcards)
    if completeness != "":
        return "--completeness %s" % (completeness)
    else:
        return ""

def get_tree_consensus_sequences(wildcards):
    pattern = get_consensus_pattern()
    consensus_sequences = [pattern.format(data_root=config["data_root"], sample=s) for s in get_sample_names()]

    if "tree_include_consensus" in config:
        consensus_sequences.append(config["tree_include_consensus"])

    return consensus_sequences

def get_tree_plot_input(wildcards):

    input_list = ["qc_analysis/{prefix}_tree.nwk", "qc_analysis/{prefix}_alleles.tsv"]
    if "assign_lineages" in config and config["assign_lineages"]:
        input_list.append("lineages/{prefix}_lineage_report.csv")
    return input_list

def get_qc_sequencing_plots(wildcards):
    prefix = get_run_name()
    out = [ "plots/%s_amplicon_covered_fraction.pdf" % (prefix),
            "plots/%s_depth_by_position.pdf" % (prefix) ]

    # add plots that need Ct values
    if get_metadata_file(wildcards) != "":
        out.append("plots/%s_amplicon_depth_by_ct.pdf" % (prefix))
    return out

def get_qc_analysis_plots(wildcards):
    prefix = get_run_name()
    out = [ "plots/%s_tree_snps.pdf" % (prefix),
            "plots/%s_amplicon_coverage_heatmap.pdf" % (prefix) ]
    return out

def get_qc_summary(wildcards):
    prefix = get_run_name()
    return "qc_reports/%s_summary_qc.tsv" % (prefix)

def get_qc_summary_metadata_opt(wildcards):
    metadata = get_metadata_file(wildcards)
    if metadata != "":
        return "--meta %s" % (metadata)
    else:
        return ""

def get_variants(wildcards):
    pattern = get_variants_pattern()
    return pattern.format(data_root=config['data_root'], sample=wildcards.sample)

def get_platform_opt(wildcards):
    if "platform" in config:
        return "--platform %s" % (config['platform'])
    else:
        return ""

def get_run_name_opt(wildcards):
    return get_run_name()

def get_consensus(wildcards):
    pattern = get_consensus_pattern()
    return pattern.format(data_root=config['data_root'], sample=wildcards.sample)

def get_run_alleles(wildcards):
    rn = get_run_name()
    return "qc_analysis/%s_alleles.tsv" % (rn)

def get_negative_control_report(wildcards):
    rn = get_run_name()
    return "qc_reports/%s_negative_control_report.tsv" % (rn)

def get_mixture_report(wildcards):
    rn = get_run_name()
    return "qc_reports/%s_mixture_report.tsv" % (rn)

def get_ambiguous_report(wildcards):
    rn = get_run_name()
    return "qc_reports/%s_ambiguous_position_report.tsv" % (rn)

def get_qc_reports(wildcards):
    out = [ get_qc_summary(wildcards) ]

    # currently these reports are only generated for illumina data
    if config.get("platform") == "illumina":
        out.append(get_mixture_report(wildcards))
        out.append(get_ambiguous_report(wildcards))

    # only try to make negative control report if NC samples have been defined
    if len(get_negative_control_samples()) > 0:
        out.append(get_negative_control_report(wildcards))
    return out

def get_report_tex_input(wildcards):
    out = get_qc_reports(wildcards)
    out.append("plots/%s_tree_snps.pdf" % (wildcards.prefix))
    if len(get_negative_control_samples()) > 0:
        out.append("plots/%s_depth_by_position_negative_control.pdf" % (wildcards.prefix))
    return out

def get_final_pdf_report(wildcards):
    rn = get_run_name()
    return "qc_reports/%s.pdf" % (rn)

def get_annotated_variants(wildcards):
    pattern = "qc_annotation/{sample}.NC_045512v2_multianno.txt"
    out = [pattern.format(sample=s) for s in get_sample_names()]
    return out

# generate the amplicon-level bed file from the input primer bed
rule make_amplicon_bed:
    input:
        primers=get_primer_bed
    output:
        "bed/amplicon.bed"
    params:
        script="primers_to_amplicons.py",
        offset=get_primer_offset,
        bed_type_opt=get_primer_bed_type_opt
    shell:
        "{params.script} --primers {input.primers} --offset {params.offset} --bed_type {params.bed_type_opt} --output {output}"

# make a bed file for the entire reference genome as a single record
# from: https://bioinformatics.stackexchange.com/questions/91/how-to-convert-fasta-to-bed
rule make_genome_bed:
    input:
        get_reference_genome_fai
    output:
        "bed/genome.bed"
    shell:
        "cat {input} | awk '{{ print $1 \"\t0\t\" $2 }}' > {output}"
