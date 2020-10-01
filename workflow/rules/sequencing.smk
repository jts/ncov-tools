#
# Rules for producing sequencing-level QC (coverage plots, genome completeness, etc)
#

# Top level rule for generating all plots
rule all_qc_sequencing:
    input:
        get_qc_sequencing_plots

# generate a bed file containing the coverage across each amplicon for a sample
rule make_amplicon_coverage:
    input:
        bam=get_bam_for_sample,
        #bed=get_amplicon_bed
        bed="bed/amplicon.bed"
    output:
        "qc_sequencing/{sample}.amplicon_coverage.bed"
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "echo -e \"reference_name\tstart\tend\tamplicon_id\tpool\tstrand\tread_count\tcovered_bases\tamplicon_length\tfraction_covered\" > {output};"
        "bedtools coverage -a {input.bed} -b {input.bam} >> {output}"

# generate a bed file containing the mean depth across each amplicon for a sample
rule make_amplicon_depth:
    input:
        bam=get_bam_for_sample,
        #bed=get_amplicon_bed
        bed="bed/amplicon.bed"
    output:
        "qc_sequencing/{sample}.amplicon_depth.bed"
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "echo -e \"reference_name\tstart\tend\tamplicon_id\tpool\tstrand\tmean_depth\" > {output};"
        "bedtools coverage -mean -a {input.bed} -b {input.bam} >> {output}"

# generate a bed file containing per-base coverage across each amplicon for a sample
rule make_amplicon_base_coverage:
    input:
        bam=get_bam_for_sample,
        #bed=get_amplicon_bed
        bed="bed/amplicon.bed"
    output:
        "qc_sequencing/{sample}.amplicon_base_coverage.bed"
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "echo -e \"reference_name\tstart\tend\tamplicon_id\tpool\tstrand\tposition\tdepth\" > {output};"
        "bedtools coverage -d -a {input.bed} -b {input.bam} >> {output}"

# generate a bed file containing per-base coverage across the entire genome
rule make_genome_per_base_coverage:
    input:
        bam=get_bam_for_sample,
        bed="bed/genome.bed"
    output:
        "{directory}/{sample}.per_base_coverage.bed"
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "echo -e \"reference_name\tstart\tend\tposition\tdepth\" > {output};"
        "bedtools coverage -d -a {input.bed} -b {input.bam} >> {output}"

# format_pileup.py requires an indexed bam file, but the connor
# pipeline doesn't create it by default. We can't assume
# write access in the data_root directory so we work around
# it by creating a directory of symlinks in the working directory
# that we can write to
rule make_tmp_bam:
    input:
        get_primer_trimmed_bam_for_sample
    output:
        "tmp_bam/{sample}.bam"
    shell:
        "ln -s \"$(readlink -f {input})\" {output}"

# build a bam index
rule make_bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    shell:
        "samtools index {input}"

# summarize samtools mpileup output in an easy to parse format
rule make_formatted_pileup:
    input:
        bam="tmp_bam/{sample}.bam",
        bai="tmp_bam/{sample}.bam.bai",
        reference=get_reference_genome
    output:
        "qc_sequencing/{sample}.fpileup.tsv"
    params:
        pileup_script = srcdir("../scripts/format_pileup.py")
    shell:
        "python {params.pileup_script} --bam {input.bam} --reference {input.reference} > {output}"

# make a file-of-filenames with the fpileup files for each sample
rule make_fpileups_fofn:
    input:
        expand("qc_sequencing/{s}.fpileup.tsv", s=get_sample_names())
    output:
        "qc_sequencing/{prefix}_fpileups.fofn",
    shell:
        'echo {input} | tr " " "\\n" > {output}'

#
# Plots
#

# heatmap of amplicon coverage (column) for each sample (row)
rule make_qc_plot_amplicon_coverage_heatmap:
    input:
        expand("qc_sequencing/{s}.amplicon_depth.bed", s=get_sample_names())
    output:
        plot="plots/{prefix}_amplicon_coverage_heatmap.pdf",
        table="qc_analysis/{prefix}_amplicon_coverage_table.tsv"
    params:
        plot_script = srcdir("../scripts/plot/plot_amplicon_coverage_heatmap.R")
    shell:
        "Rscript {params.plot_script} --path qc_sequencing --output {output.plot} --table {output.table}"

def get_metadata_opt(wildcards):
    metadata_file = get_metadata_file(wildcards)
    if metadata_file != "":
        return "-m %s" % get_metadata_file(wildcards)
    else:
        return ""

# plot coverage along the genome for each sample
rule make_qc_plot_depth_by_position:
    input:
        expand("qc_sequencing/{s}.per_base_coverage.bed", s=get_sample_names())
    output:
        "plots/{prefix}_depth_by_position.pdf"
    params:
        plot_script = srcdir("../scripts/plot/plot_qc_sequencing.R"),
        metadata_opt=get_metadata_opt
    shell:
        "Rscript {params.plot_script} -t depth_by_position -o {output} {params.metadata_opt}"

# plot coverage, except only for the negative controls
rule make_qc_plot_depth_by_position_negative_controls:
    input:
        expand("qc_sequencing_negative_control/{s}.per_base_coverage.bed", s=get_negative_control_samples())
    output:
        "plots/{prefix}_depth_by_position_negative_control.pdf"
    params:
        plot_script = srcdir("../scripts/plot/plot_qc_sequencing.R"),
        metadata_opt=get_metadata_opt
    shell:
        "Rscript {params.plot_script} -t negative_control_depth_by_position -o {output} {params.metadata_opt}"

# 
rule make_qc_plot_amplicon_depth_by_ct:
    input:
        files=expand("qc_sequencing/{s}.amplicon_depth.bed", s=get_sample_names()),
        metadata=get_metadata_file
    output:
        "plots/{prefix}_amplicon_depth_by_ct.pdf"
    params:
        plot_script = srcdir("../scripts/plot/plot_qc_sequencing.R")
    shell:
        "Rscript {params.plot_script} -t amplicon_depth_by_ct -o {output} -m {input.metadata}"

#
rule make_qc_plot_fraction_covered_by_amplicon:
    input:
        expand("qc_sequencing/{s}.amplicon_coverage.bed", s=get_sample_names())
    output:
        "plots/{prefix}_amplicon_covered_fraction.pdf"
    params:
        plot_script = srcdir("../scripts/plot/plot_qc_sequencing.R")
    shell:
        "Rscript {params.plot_script} -t amplicon_covered_fraction -o {output}"

#
rule make_qc_genome_completeness_by_ct:
    input:
        qc="qc_analysis/merged.qc.csv",
        metadata=get_metadata_file
    output:
        "plots/{prefix}_genome_completeness_by_ct.pdf"
    params:
        plot_script = srcdir("../scripts/plot/plot_qc_sequencing.R")
    shell:
        "Rscript {params.plot_script} -t genome_completeness_by_ct -o {output} -m {input.metadata}"
