#! /usr/bin/python

import argparse
from datetime import datetime
import glob
import csv
import sys
import os
import re

#
# utility class to assist in converting
# a tsv into a latex table - it can
# rename header columns (name_map)
# perform arbitrary transforms (row_func)
# and filter out columns
#
class TableFormatter:
    def __init__(self):
        self.name_map = dict()
        self.row_func = dict()
        self.column_filter = dict()
        self.table_spec = ""
        self.size = "normalsize"

# convert a pdf to a collection of pngs, returning a list of the generated filenames
def pdf_to_png(pdf_name, output_directory="report_images"):

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    prefix = os.path.basename(pdf_name).replace(".pdf", "")
    os.system("pdftoppm %s %s/%s -png" % (pdf_name, output_directory, prefix))
    return sorted(glob.glob("%s/%s*.png" % (output_directory, prefix)))

# latex starting boilerplate to set up the document class, packages, etc
def write_preamble():

    p = r'''
\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage{graphicx}
\usepackage{array}
\usepackage{float}
\usepackage{longtable}
\usepackage[margin=0.15in]{geometry}
\newcolumntype{C}[1]{>{\centering\let\newline\\\arraybackslash\hspace{0pt}}m{#1}}
\begin{document}'''
    print(p)

# latex ending boilerplate
def write_postamble():

    p = r'''
\end{document}'''
    print(p)

# count the number of rows in a tsv file
def count_tsv(filename):

    c = 0
    with(open(filename)) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            c += 1
    return c

def format_ct(s):

    try:
        return "%.1f" % (float(s))
    except:
        return "NA"

def escape_latex(s):
    return s.replace("_", "\_")

# high-level function to transform a tsv into a latex table
# performing remapping of column names and arbitrary transformation
# of values within each column (e.g. escaping characters latex doesn't like)
# using table formatter
def tsv_to_table(filename, table_formatter):
    with(open(filename)) as f:
        reader = csv.DictReader(f, delimiter="\t")

        rows = list()
        header = list()
        
        for row in reader:

            # remove columns
            for k in table_formatter.column_filter:
                del row[k]

            if len(header) == 0:

                # remap column names
                for k in row.keys():
                    if k in table_formatter.name_map:
                        header.append(table_formatter.name_map[k])
                    else:
                        header.append(escape_latex(k))
                
            # transform with row func
            for k in row:
                if k in table_formatter.row_func:
                    row[k] = table_formatter.row_func[k](row[k])
            rows.append(row.values())

        # write latex to stdout
        write_table(table_formatter.table_spec, header, rows, table_formatter.size)

# latex for displaying an image
def write_image(image_filename, scale=1.0):
    print(r"\begin{center}")
    print("\includegraphics[scale=%f]{%s}" % (scale, image_filename))
    print(r"\end{center}")

# latex for displaying a table 
def write_table(spec, header, rows, size):
    print(r"\begin{center}")
    print(r"\%s" % size)

    print(r"\begin{longtable}{%s}" % (spec))
    print(r"\hline")
    print(" & ".join(header) + r" \\ \hline")
    print(r"\endhead")
    for r in rows:
        print(" & ".join([escape_latex(v) for v in r]) + r" \\ \hline")

    print(r"\end{longtable}")
    print(r"\normalsize")
    print(r"\end{center}")

# parse a sample name from a filename, assuming the
# sample name is the first entry
def filename_to_sample(fn):
    return os.path.basename(fn).split(".")[0]

# write the section of the report describing the negative control checks
def write_negative_control_section():

    print(r"\section{Negative Control}")


    plot_path = args.negative_control_depth_figure.format(run_name=args.run_name)

    if os.path.exists(plot_path):
        nc_depth_figures = pdf_to_png(plot_path)

        # plot each depth image for the negative control samples
        for depth_png in nc_depth_figures:
            write_image(depth_png, 0.75)
    else:
        print(r"Negative control depth plots not found.")


    # write summary table
    table_path = args.negative_control_table.format(run_name=args.run_name)
    
    if os.path.exists(table_path):

        # set up TableFormatter to transform the tsv into a nicer display
        tf = TableFormatter()

        # map to rename header columns into more interpretable names
        tf.name_map = { "file" : "Sample",
                        "genome_covered_bases" : "Covered",
                        "genome_total_bases" : "Target Footprint",
                        "genome_covered_fraction" : "Percent Covered",
                        "amplicons_detected" : "Amplicons Detected" 
                      }

        # map to transform selected row columns into more readable values
        tf.row_func = { "file" : lambda value : filename_to_sample(value),
                        "amplicons_detected" : lambda value : value.replace(",", ", "), 
                        "genome_covered_fraction" : lambda value : "%.1f" % (float(value) * 100.0)
                      }

        tf.table_spec = "|c|c|c|c|c|p{5cm}|"
        tsv_to_table(table_path, tf)
    else:
        print(r"Negative control summary table not found.")
    return

# write the section of the report containing the tree-snps plot
# and the ambiguity report
def write_tree_section():
    print("\section{Sequence Variation}")
   
    tree_plot_path = args.tree_figure.format(run_name=args.run_name)
    write_image(tree_plot_path, 0.4)

    # Ambiguity subsection
    print(r"\subsection{Ambiguity Report}")

    if args.platform == "illumina":
        ambiguity_filename = args.ambiguity_table.format(run_name=args.run_name)
        row_count = count_tsv(ambiguity_filename)

        if row_count == 0:
            print(r"No positions in the genome had an ambiguous consensus base (IUPAC, but not N) in multiple samples.")
        else:
            s = "This table reports positions in the genome that had an ambiguous consensus base"\
                " (IUPAC, but not N) in multiple samples. This can be evidence of contamination"\
                " so these positions should be investigated."
            print(s)

            # Make the ambiguity table
            tf = TableFormatter()
            tf.name_map = { "position" : "Genome Position",
                            "count"    : "Number of ambiguous samples",
                            "alleles"  : "Observed Alleles" }
            tf.row_func = dict()
            tf.table_spec = "{|c|c|c|}"
            tsv_to_table(args.ambiguity_table.format(run_name=args.run_name), tf)
    else:
        print(r"The ambiguity report is not currently supported for Oxford Nanopore data")

    # Mixture/Contamination subsection
    print(r"\subsection{Mixture Report}")
    
    # The mixture table can be quite large so we report the samples as a list
    if args.platform == "illumina":
        mixture_report_fn = args.mixture_table.format(run_name=args.run_name)
        mixture_samples = set()
        with(open(mixture_report_fn)) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                mixture_samples.add(row['sample_a'])
        
        if len(mixture_samples) > 0:
            s = "The following samples were detected by \\texttt{mixture\_report.py} as having"\
                " read evidence for multiple distinct sequences. These samples should be checked"\
                " for contamation: "
            print(s + ", ".join([ "\\textbf{%s}" % escape_latex(a) for a in mixture_samples]) + ".")
        else:
            print(r"No samples were detected by \texttt{mixture\_report.py} as having read evidence for multiple distinct sequences.")
    else:
        print("The mixture report is not currently supported for Oxford Nanopore data")

# write the large per-sample QC table
def write_summary_qc_section():
    print(r"\section{Sample-level QC}")
    print(r"This table contains QC metrics and warning flags for each sample within this sequencing run.")

    tf = TableFormatter
    tf.size = "scriptsize"
    tf.name_map = { "sample" : "Sample",
                    "num_consensus_snvs" :  "Consensus SNVs",
                    "num_consensus_iupac" : "Consensus IUPAC",
                    "num_variants_snvs" : "Variant SNVs",
                    "num_variants_indel" : "Variant Indels",
                    "num_variants_indel_triplet" : "Variant Triplet Indels",
                    "qpcr_ct" : "ct",
                    "collection_date" : "Date",
                    "genome_completeness" : "Percent Complete",
                    "qc_pass" : "QC flags" }

    tf.row_func = { "qc_pass" : lambda value : value.replace(",", ", ").replace("POSSIBLE_FRAMESHIFT_INDELS", "POSSIBLE_FRAMESHIFT"),
                    "genome_completeness" : lambda value : "%.1f" % (float(value) * 100.0),
                    "qpcr_ct" : lambda value : format_ct(value) }

    tf.column_filter = [ "run_name", "mean_sequencing_depth", 
                         "median_sequencing_depth", "num_consensus_n", 
                         "num_weeks", "scaled_variants_snvs" ]
    tf.table_spec = "{|c|C{1.3cm}|C{1.3cm}|C{1.0cm}|C{1.0cm}|C{1.0cm}|c|c|C{1.2cm}|C{4.0cm}|}"
    tsv_to_table(args.summary_qc_table.format(run_name=args.run_name), tf)

# parse arguments, defaults are based on ncov-tools paths/filenames
parser = argparse.ArgumentParser()
parser.add_argument('--negative-control-depth-figure', type=str, default="plots/{run_name}_depth_by_position_negative_control.pdf")
parser.add_argument('--negative-control-table', type=str, default="qc_reports/{run_name}_negative_control_report.tsv")
parser.add_argument('--tree-figure', type=str, default="plots/{run_name}_tree_snps.pdf")
parser.add_argument('--ambiguity-table', type=str, default="qc_reports/{run_name}_ambiguous_position_report.tsv")
parser.add_argument('--mixture-table', type=str, default="qc_reports/{run_name}_mixture_report.tsv")
parser.add_argument('--summary-qc-table', type=str, default="qc_reports/{run_name}_summary_qc.tsv")
parser.add_argument('--negative-control-tsv', type=str, default="")
parser.add_argument('--run-name', type=str, default="", required=True)
parser.add_argument('--platform', type=str, default="", required=True)
args = parser.parse_args()

#
# Generate report
#
write_preamble()

print(r"\section{Summary}")
print("QC report generated for %s on %s by ncov-tools" % (escape_latex(args.run_name), datetime.today().strftime('%Y-%m-%d')))

write_tree_section()

write_negative_control_section()

write_summary_qc_section()

write_postamble()


