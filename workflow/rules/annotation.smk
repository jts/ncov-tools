#
# Helper functions
#
def get_aa_table(wildcards):
    pattern = "qc_annotation/{sample}_aa_table.tsv"
    out = [pattern.format(sample=s) for s in get_sample_names()]
    return out

def get_vcf_file(wildcards):
    if config['platform'] == 'illumina':
        pattern = "qc_annotation/{sample}.pass.vcf.gz"
    elif config['platform'] == 'oxford-nanopore':
        pattern = "data/{sample}.pass.vcf.gz"
    return pattern

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

def get_java_heap_opt():
    '''
    Allocation heap space to the JVM for use with SNPEff
    '''
    return config.get('snpeff_java_heap', '8g')

#
# Rules for annotating variants with functional consequence
#
rule annotate_snpeff:
    input:
        get_recurrent_heatmap_plot

rule convert_ivar_to_vcf:
    input:
        "data/{sample}.variants.tsv"
    output:
        "qc_annotation/{sample}.pass.vcf"
    params:
        script=srcdir("../scripts/ivar_variants_to_vcf.py")
    shell:
        "{params.script} {input} {output}"

rule compress_vcf:
    input:
        "qc_annotation/{sample}.pass.vcf"
    output:
        "qc_annotation/{sample}.pass.vcf.gz"
    shell:
        "gzip {input}" 

rule run_snpeff:
    input:
        get_vcf_file
    output:
        "qc_annotation/{sample}.ann.vcf"
    params:
        jar=os.environ["SNPEFFJAR"],
        db="MN908947.3",
        java_heap=get_java_heap_opt
    shell:
        "java -Xmx{params.java_heap} -jar {params.jar} {params.db} {input} > {output}"


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

