import os
import re
#
# Helper functions
#

def get_aa_table(wildcards):
    pattern = "qc_annotation/{sample}_aa_table.tsv"
    out = [pattern.format(sample=s) for s in get_sample_names()]
    return out

def get_vcf_file(wildcards):
    data_root = config['data_root']
    variants_pattern = get_variants_pattern()
    is_vcf = variants_pattern.endswith(".vcf") or variants_pattern.endswith(".vcf.gz")
    
    if is_vcf:
        pattern = variants_pattern.format(data_root=config['data_root'], sample=wildcards.sample)
    else:
        # assume illumina ivar pipeline, this will trigger conversion
        assert(config['platform'] == 'illumina')
        pattern = "qc_annotation/{sample}.pass.vcf.gz"
    return pattern

def get_variants_file(wildcards):
    pattern = get_variants_pattern()
    out = pattern.format(data_root=config['data_root'], sample=wildcards.sample)
    return out

def get_snpeff_vcf_files(wildcards):
    pattern = "qc_annotation/{sample}.ann.vcf"
    out = [pattern.format(sample=s) for s in get_sample_names()]
    return out

def get_threshold_opt(wildcards):
    threshold = config.get("rec_threshold", "2")
    return threshold

def get_recurrent_heatmap_plot(wildcards):
    prefix = get_run_name()
    out = "plots/%s_aa_mutation_heatmap.pdf" % (prefix)
    return out

def get_primer_snp_depth_table(wildcards):
    prefix = get_run_name()
    out = "qc_annotation/%s.primer_snp_depth_summary.tsv" % (prefix)
    return out

#
# Rules for annotating variants with functional consequence
#
rule all_qc_annotation:
    input:
        get_recurrent_heatmap_plot,
        get_primer_snp_depth_table

rule build_snpeff_db:
    input:
        expand(os.environ['CONDA_PREFIX'] + '/share/snpeff-5.0-0/data/MN908947.3/snpEffectPredictor.bin')

rule download_db_files:
    output:
        expand(os.environ['CONDA_PREFIX'] + '/share/snpeff-5.0-0/data/MN908947.3/snpEffectPredictor.bin')
    params:
        script=srcdir("../scripts/build_db.py")
    shell:
        "python {params.script}"

rule convert_ivar_to_vcf:
    input:
        get_variants_file
    output:
        "qc_annotation/{sample}.pass.vcf"
    params:
        script=srcdir("../scripts/ivar_variants_to_vcf.py")
    shell:
        "{params.script} {input} {output}"

rule compress_vcf:
    input:
        "qc_annotation/{prefix}.vcf"
    output:
        "qc_annotation/{prefix}.vcf.gz"
    shell:
        "gzip {input}" 

rule run_snpeff:
    input:
        get_vcf_file
    output:
        "qc_annotation/{sample}.ann.vcf"
    params:
        script="snpEff",
        db="MN908947.3",
        aa_letter="-hgvs1LetterAa"
    shell:
        "{params.script} {params.aa_letter} {params.db} {input} > {output}"


rule convert_annotated_vcf_to_aa_table:
    input:
        vcf="qc_annotation/{sample}.ann.vcf"
    output:
        "qc_annotation/{sample}_aa_table.tsv"
    params:
        script=srcdir("../scripts/convert_recurrent_snpeff_data.py"),
        sample="{sample}"
    shell:
        "python {params.script} --file {input.vcf} --sample {params.sample} --output {output}"

rule create_recurrent_mutation_heatmap:
    input:
        get_aa_table
    output:
        "plots/{prefix}_aa_mutation_heatmap.pdf"
    params:
        script=srcdir("../scripts/plot/plot_recurrent_variant_heatmap_snpeff.R"),
        threshold=get_threshold_opt
    shell:
        "Rscript {params.script} --path qc_annotation --output {output} --threshold {params.threshold}"

rule create_primer_snp_bed:
    input:
        vcf=get_vcf_file,
        primer_bed=get_primer_bed
    output:
        "qc_annotation/{sample}.primer_snps.bed"
    shell:
        "bedtools intersect -a {input.vcf} -b {input.primer_bed} -wb > {output}"

rule create_primer_snp_depth:
    input:
        primer_snps="qc_annotation/{sample}.primer_snps.bed",
        amplicon_depth="qc_sequencing/{sample}.amplicon_depth.bed"
    output:
        "qc_annotation/{sample}.primer_snp_depth.tsv"
    params:
        script=srcdir("../scripts/primer_snp_depth.py")
    shell:
        "python {params.script} --sample-name {wildcards.sample} --primer-snp-bed {input.primer_snps} --amplicon-depth-tsv {input.amplicon_depth} > {output}"

rule merge_primer_snp_depth:
    input:
        expand("qc_annotation/{s}.primer_snp_depth.tsv", s=get_sample_names())
    output:
        "qc_annotation/{prefix}.primer_snp_depth_merged.tsv"
    shell:
        "cat {input} | awk 'NR == 1 || $0 !~ /primer_name/' > {output}"

rule summarize_primer_snp_depth:
    input:
        "qc_annotation/{prefix}.primer_snp_depth_merged.tsv"
    output:
        "qc_annotation/{prefix}.primer_snp_depth_summary.tsv"
    params:
        script=srcdir("../scripts/primer_snp_summary.py")
    shell:
        "python {params.script} --input-tsv {input} > {output}"

