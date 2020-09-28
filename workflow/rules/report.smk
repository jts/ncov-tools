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

rule make_report_tex:
    input:
        get_report_tex_input
    output:
        "qc_reports/{prefix}.tex"
    params:
        script=srcdir("../scripts/generate_report.py"),
        run_name=get_run_name,
        platform_opt=get_platform_opt
    shell:
        "python {params.script} --run-name {params.run_name} {params.platform_opt} > {output}"

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
