def get_threshold_opt(wildcards):
    if "recurrent_threshold" in config:
        return f'--threshold {config["recurrent_threshold"]}'
    else:
        return ''

def get_aa_heatmap_basesize_opt(wildcards):
    if 'basesize_recurrent_heatmap' in config:
        return f'--base_size {config["basesize_recurrent_heatmap"]}'
    else:
        return ''

def get_aa_heatmap(wildcards):
    prefix = get_run_name()
    return f'plots/{prefix}_aa_mutation_heatmap.pdf'

def get_annotated_variants(wildcards):
    pattern = "qc_annotation/{sample}.NC_045512v2_multianno.txt"
    out = [pattern.format(sample=s) for s in get_sample_names()]
    return out

def get_sarscov2db_opt(wildcards):
    return config.get('sarscov2db', '')

#
# Rules for annotating variants with functional consequence
#
rule annotate_variants:
    input:
        get_aa_heatmap

# run convert a variants file into annovar's input format
rule make_annovar_input:
    input:
        samplevariant=get_variants
    output:
        "qc_annotation/{sample}.avinput"
    params:
        script=srcdir("../scripts/create_annovar_input.py")
    shell:
        "python {params.script} --file {input.samplevariant} --output {output}"

# run annovar to predict the functional consequence of each mutation
rule run_table_annovar:
    input:
        "qc_annotation/{sample}.avinput"
    output:
        "qc_annotation/{sample}.NC_045512v2_multianno.txt"
    params:
        script="table_annovar.pl",
        buildver="NC_045512v2",
        sarscov2db=get_sarscov2db_opt,
        outfile="qc_annotation/{sample}"
    shell:
        "{params.script} --buildver {params.buildver} {input} {params.sarscov2db} -protocol avGene -operation g --remove --otherinfo --outfile {params.outfile}"

rule make_aa_heatmap:
    input:
        get_annotated_variants
    output:
        "plots/{prefix}_aa_mutation_heatmap.pdf"
    params:
        script=srcdir("../scripts/plot/plot_recurrent_variant_heatmap.R"),
        threshold_opt=get_threshold_opt,
        base_size_opt=get_aa_heatmap_basesize_opt
    shell:
        "Rscript {params.script} --path qc_annotation/ --output {output} {params.threshold_opt} {params.base_size_opt}"

