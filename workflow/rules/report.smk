#
# Rules for producing cancogen QC reports (negative control, contamination checks, etc) 
#
rule all_qc_summary:
    input:
        get_qc_summary

rule all_qc_negative_control:
    input:
        get_negative_control_report

rule all_mixture_report:
    input:
        get_mixture_report

rule all_ambiguous_report:
    input:
        get_ambiguous_report

rule all_qc_reports:
    input:
        get_qc_reports

rule all_final_report:
    input:
        get_final_pdf_report

# generate the negative control report, containing coverage statistics
# for each negative control sample in the config
rule make_negative_control_report:
    input:
        bed=get_negative_control_bed
    output:
        "qc_reports/{prefix}_negative_control_report.tsv"
    params:
        script=srcdir("../scripts/negative_control_check.py"),
        primer_prefix=get_primer_prefix
    shell:
        "python {params.script} --primer_prefix {params.primer_prefix} {input.bed} > {output}"

# detect possible contaimination or mixtures of samples based on allele frequencies
rule make_mixture_report:
    input:
        fpileups="qc_sequencing/{prefix}_fpileups.fofn",
        alleles="qc_analysis/{prefix}_alleles.tsv"
    output:
        "qc_reports/{prefix}_mixture_report.tsv"
    params:
        script=srcdir("../scripts/mixture_check.py")
    shell:
        "python {params.script} --fpileup {input.fpileups} --alleles {input.alleles} > {output}"

# make a simple report containing positions that frequently contain IUPAC ambiguity codes
rule make_ambiguous_position_report:
    input:
        alleles="qc_analysis/{prefix}_alleles.tsv"
    output:
        "qc_reports/{prefix}_ambiguous_position_report.tsv"
    params:
        script=srcdir("../scripts/ambiguous_position_check.py")
    shell:
        "python {params.script} --alleles {input.alleles} --min-count 3 > {output}"

# generate the .tex for the single PDF report containing all QC reports
rule make_report_tex:
    input:
        get_report_tex_input
    output:
        "qc_reports/{prefix}.tex"
    params:
        script=srcdir("../scripts/generate_report.py"),
        run_name=get_run_name,
        platform_opt=get_platform_opt,
        voc_lineages=get_voc_pango_lineages
    shell:
        "python {params.script} --run-name {params.run_name} --voc-lineages {params.voc_lineages} {params.platform_opt} > {output}"

# generate the PDF from the .tex report
rule make_report_pdf:
    input:
        "qc_reports/{prefix}.tex"
    output:
        "qc_reports/{prefix}.pdf"
    shell:
        # The long table package needs pdflatex to be run up to 3 times
        # to correctly typeset the tables
        "pdflatex -output-directory qc_reports {input} &&"
        "pdflatex -output-directory qc_reports {input} &&"
        "pdflatex -output-directory qc_reports {input}"
