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

def get_snpeff_dirs():
    snpeff_dirs = list()
    for snpeff_dir in os.scandir('/'.join([os.environ['CONDA_PREFIX'], 'share'])):
        if snpeff_dir.name.startswith('snpeff'):
            snpeff_dirs.append(snpeff_dir.name)
    return snpeff_dirs

def get_bcftools_vcf_files(wildcards):
    """
    Return a list of bcftools csq VCF files
    """
    pattern = "qc_bcftools/{sample}.csq.vcf.gz"
    out = [pattern.format(sample=s) for s in get_sample_names()]
    return out

def get_bcftools_aa_table_files(wildcards):
    """
    Returns a list of aa table files from bcftools csq VCF files
    """
    pattern = "qc_bcftools/{sample}_aa_table.tsv"
    out = [pattern.format(sample=s) for s in get_sample_names()]
    return out

def get_bcftools_recurrent_heatmap_plot(wildcards):
    prefix = get_run_name()
    out = "plots/%s_aa_mutation_heatmap_bcftools.pdf" % (prefix)
    return out

#
# Rules for annotating variants with functional consequence
#
rule run_bcftools:
    input:
        get_bcftools_recurrent_heatmap_plot
        #get_bcftools_aa_table_files
        #get_bcftools_vcf_files

rule all_qc_annotation:
    input:
        get_recurrent_heatmap_plot

rule build_snpeff_db:
    input:
        expand(os.environ['CONDA_PREFIX'] + '/share/{snpeff_dir}/data/MN908947.3/snpEffectPredictor.bin', snpeff_dir=get_snpeff_dirs())

rule download_db_files:
    input:
        expand(os.environ['CONDA_PREFIX'] + '/share/{snpeff_dir}', snpeff_dir=get_snpeff_dirs())
    output:
        expand(os.environ['CONDA_PREFIX'] + '/share/{snpeff_dir}/data/MN908947.3/snpEffectPredictor.bin', snpeff_dir=get_snpeff_dirs())
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
        "bgzip {input} && tabix -p vcf {output}" 

rule run_snpeff:
    input:
        get_vcf_file
    output:
        "qc_annotation/{sample}.ann.vcf"
    params:
        script="snpEff",
        db="MN908947.3",
        aa_letter="-hgvs1LetterAa",
        no_log="-noLog"
    shell:
        "{params.script} {params.no_log} {params.aa_letter} {params.db} {input} > {output}"


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


rule run_bcftools_csq:
    input:
        get_vcf_file
    output:
        "qc_bcftools/{sample}.csq.vcf.gz"
    params:
        script="bcftools csq",
        ref=get_reference_genome,
        gff=get_annotation_gff,
        outtype="z"
    shell:
        "{params.script} --fasta-ref {params.ref} --gff-annot {params.gff} --output {output} --output-type {params.outtype} {input}"

rule convert_bcftools_vcf_to_aa_table:
    input:
        vcf="qc_bcftools/{sample}.csq.vcf.gz"
    output:
        "qc_bcftools/{sample}_aa_table.tsv"
    params:
        script=srcdir("../scripts/convert_recurrent_bcftools_data.py"),
        sample="{sample}"
    shell:
        "python {params.script} --vcf {input.vcf} --sample {params.sample} --output {output}"

rule create_bcftools_recurrent_mutation_heatmap:
    input:
        get_bcftools_aa_table_files
    output:
        "plots/{prefix}_aa_mutation_heatmap_bcftools.pdf"
    params:
        script=srcdir("../scripts/plot/plot_recurrent_variant_heatmap_bcftools.R"),
        threshold=get_threshold_opt
    shell:
        "Rscript {params.script} --path qc_bcftools --output {output} --threshold {params.threshold}"
