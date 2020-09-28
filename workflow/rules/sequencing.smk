rule all_qc_sequencing:
    input:
        get_qc_sequencing_plots

#
# generate coverage QC data using bedtools
#
rule make_amplicon_mean_coverage:
    input:
        bam=get_bam_for_sample,
        #bed=get_amplicon_bed
        bed="bed/amplicon.bed"
    output:
        "qc_sequencing/{sample}.mean_coverage.bed"
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "echo -e \"reference_name\tstart\tend\tamplicon_id\tpool\tstrand\tmean_coverage\" > {output};"
        "bedtools coverage -mean -a {input.bed} -b {input.bam} >> {output}"

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

# https://bioinformatics.stackexchange.com/questions/91/how-to-convert-fasta-to-bed
rule make_genome_bed:
    input:
        get_reference_genome_fai
    output:
        "bed/genome.bed"
    shell:
        "cat {input} | awk '{{ print $1 \"\t0\t\" $2 }}' > {output}"

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

# pysam's index_filename option is broken so we have
# to do some hacky symlinking to work around it
rule make_tmp_bam:
    input:
        get_primer_trimmed_bam_for_sample
    output:
        "tmp_bam/{sample}.bam"
    shell:
        "ln -s \"$(readlink -f {input})\" {output}"

rule make_bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    shell:
        "samtools index {input}"

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

rule make_negative_control_report:
    input:
        bed=get_negative_control_bed
    output:
        "qc_reports/{prefix}_negative_control_report.tsv"
    params:
        script=srcdir("../scripts/negative_control_check.py")
    shell:
        "python {params.script} {input.bed} > {output}"

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

rule make_qc_plot_depth_by_position:
    input:
        expand("qc_sequencing/{s}.per_base_coverage.bed", s=get_sample_names())
    output:
        "plots/{prefix}_depth_by_position.pdf"
    params:
        plot_script = srcdir("../scripts/plot/plot_qc_sequencing.R"),
        metadata=get_metadata_file
    shell:
        "Rscript {params.plot_script} depth_by_position {wildcards.prefix} {params.metadata}"

rule make_qc_plot_depth_by_position_negative_controls:
    input:
        expand("qc_sequencing_negative_control/{s}.per_base_coverage.bed", s=get_negative_control_samples())
    output:
        "plots/{prefix}_depth_by_position_negative_control.pdf"
    params:
        plot_script = srcdir("../scripts/plot/plot_qc_sequencing.R"),
        metadata=get_metadata_file
    shell:
        "Rscript {params.plot_script} negative_control_depth_by_position {wildcards.prefix} {params.metadata}"

rule make_qc_plot_amplicon_depth_by_ct:
    input:
        files=expand("qc_sequencing/{s}.amplicon_depth.bed", s=get_sample_names()),
        metadata=get_metadata_file
    output:
        "plots/{prefix}_amplicon_depth_by_ct.pdf"
    params:
        plot_script = srcdir("../scripts/plot/plot_qc_sequencing.R")
    shell:
        "Rscript {params.plot_script} amplicon_depth_by_ct {wildcards.prefix} {input.metadata}"

rule make_qc_plot_fraction_covered_by_amplicon:
    input:
        expand("qc_sequencing/{s}.amplicon_coverage.bed", s=get_sample_names())
    output:
        "plots/{prefix}_amplicon_covered_fraction.pdf"
    params:
        plot_script = srcdir("../scripts/plot/plot_qc_sequencing.R")
    shell:
        "Rscript {params.plot_script} amplicon_covered_fraction {wildcards.prefix}"

rule make_qc_genome_completeness_by_ct:
    input:
        qc="qc_analysis/merged.qc.csv",
        metadata=get_metadata_file
    output:
        "plots/{prefix}_genome_completeness_by_ct.pdf"
    params:
        plot_script = srcdir("../scripts/plot/plot_qc_sequencing.R")
    shell:
        "Rscript {params.plot_script} genome_completeness_by_ct {wildcards.prefix} {input.metadata}"
