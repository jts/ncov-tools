#
# Rules for annotating variants with functional consequence
#
rule annotate_variants:
    input:
        get_annotated_variants

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
