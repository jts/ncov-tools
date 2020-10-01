#
# Rules for producing analysis-level QC (phylogenetic trees, statistics on number of called variants)
#

# top-level rule for generating all plots 
rule all_qc_analysis:
    input:
        get_qc_analysis_plots

#
# Make a tree
#

# merge the consensus sequences together into one fasta, applying a
# completeness threshold to avoid including very poor samples
rule make_merged_consensus:
    input:
        get_tree_consensus_sequences
    output:
        "qc_analysis/{prefix}_consensus.fasta"
    params:
        completeness_opt=get_completeness_threshold_opt,
        rename_script = srcdir("../scripts/preprocess_consensus.py")
    shell:
        "python {params.rename_script} {params.completeness_opt} {input} > {output}"

# make the multiple sequence alignment
rule make_msa:
    input:
        consensus="qc_analysis/{prefix}_consensus.fasta",
        reference=get_reference_genome
    output:
        "qc_analysis/{prefix}_aligned.fasta"
    shell:
        "augur align --sequences {input.consensus} --reference-sequence {input.reference} --output {output} --fill-gaps"

# make the basic tree from the MSA
rule make_tree_raw:
    input:
        "qc_analysis/{prefix}_aligned.fasta"
    output:
        "qc_analysis/{prefix}_tree_raw.nwk"
    shell:
        "augur tree --alignment {input} --output {output}"

# reroot the tree to the reference genome
rule make_tree_final:
    input:
        tree="qc_analysis/{prefix}_tree_raw.nwk",
        reference=get_reference_genome
    output:
        "qc_analysis/{prefix}_tree.nwk"
    shell:
        "nw_reroot {input} `head -1 {input.reference} | tr -d \">\"` > {output}"

# parse the MSA to identify differences wrt the reference
rule make_alleles:
    input:
        "qc_analysis/{prefix}_aligned.fasta"
    output:
        "qc_analysis/{prefix}_alleles.tsv"
    params:
        alleles_script = srcdir("../scripts/align2alleles.py")
    shell:
        "python {params.alleles_script} --reference-name MN908947.3 {input} > {output}"

# assign lineages to the consensus genomes using pangolin
rule make_lineage_assignments:
    input:
        "qc_analysis/{prefix}_consensus.fasta"
    output:
        "lineages/{prefix}_lineage_report.csv"
    threads: workflow.cores
    shell:
        "pangolin -t {threads} --outfile {output} {input}"

# make the primary plot, containing the phylogenetic tree and associated mutations
rule make_qc_tree_snps:
    input: get_tree_plot_input
    output:
        "plots/{prefix}_tree_snps.pdf"
    params:
        plot_script = srcdir("../scripts/plot/plot_tree_snps.R")
    shell:
        "Rscript {params.plot_script} {output} {input}"

# generate the summary QC metrics for each sample
rule make_sample_qc_summary:
    input:
        alleles=get_run_alleles,
        samplecoverage="qc_sequencing/{sample}.per_base_coverage.bed",
        samplevariants=get_variants,
        sampleconsensus=get_consensus
    output:
        "qc_analysis/{sample}.summary.qc.tsv"
    params:
        py_script="get_qc.py",
        metadata_opt=get_qc_summary_metadata_opt,
        platform_opt=get_platform_opt,
        run_name_opt=get_run_name_opt
    shell:
        "{params.py_script} --alleles {input.alleles} \
                            --coverage {input.samplecoverage} \
                            --variants {input.samplevariants} {params.metadata_opt} \
                            --indel \
                            --consensus {input.sampleconsensus} {params.platform_opt} \
                            --sample {wildcards.sample} \
                            --run_name {params.run_name_opt} > {output}"

# merge the per-sample summary files into the run-level report
rule make_full_qc_summary:
    input:
        expand("qc_analysis/{s}.summary.qc.tsv", s=get_sample_names())
    output:
        "qc_reports/{prefix}_summary_qc.tsv"
    params:
        py_script="collect_qc_summary.py"
    shell:
        "{params.py_script} --path qc_analysis > {output}"

# merge the nextflow QC output csv files
rule merge_artic_qc:
    input:
        expand(config["data_root"] + "/{s}.qc.csv", s=get_sample_names())
    output:
        "qc_analysis/merged.qc.csv"
    shell:
        "cat {input} | awk 'NR == 1 || $0 !~ /qc_pass/' > {output}"
