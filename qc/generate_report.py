#! /usr/bin/python

import argparse
from datetime import datetime
import glob
import csv
import sys
import os
import re

# convert a pdf to a collection of pngs, returning a list of filenames
def pdf_to_png(pdf_name, output_directory="report_images"):

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    prefix = os.path.basename(pdf_name).replace(".pdf", "")
    os.system("pdftoppm %s %s/%s -png" % (pdf_name, output_directory, prefix))
    return sorted(glob.glob("%s/%s*.png" % (output_directory, prefix)))

# latex starting boilerplate
def write_preamble():

    p = r'''
\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage{graphicx}
\usepackage[margin=0.15in]{geometry}
\begin{document}'''
    print(p)

# latex ending boilerplate
def write_postamble():

    p = r'''
\end{document}'''
    print(p)

def escape_latex(s):
    return s.replace("_", "\_")

#
def write_image(image_filename, scale=1.0):
    print(r"\begin{center}")
    print("\includegraphics[scale=%f]{%s}" % (scale, image_filename))
    print(r"\end{center}")

def write_table(spec, header, rows):
    print(r"\begin{center}")
    print(r"\begin{tabular}{%s}" % (spec))
    print(r"\hline")
    print(" & ".join(header) + r" \\ \hline")
    for r in rows:
        print(" & ".join([escape_latex(v) for v in r]) + r" \\ \hline")

    print(r"\end{tabular}")
    print(r"\end{center}")

def filename_to_sample(fn):
    return os.path.basename(fn).split(".")[0]

def write_negative_control_section():

    plot_path = args.negative_control_depth_figure.format(run_name=args.run_name)

    nc_depth_figures = pdf_to_png(plot_path)

    print(r"\section{Negative Control}")

    # Plot each depth image for the negative control samples
    for depth_png in nc_depth_figures:
        write_image(depth_png, 0.75)

    table_path = args.negative_control_table.format(run_name=args.run_name)
    
    # Load negative control table
    with(open(table_path)) as f:
        reader = csv.DictReader(f, delimiter="\t")

        rows = list()
        header = list()

        name_map = { "genome_covered_bases" : "Covered",
                     "genome_total_bases" : "Target Footprint",
                     "genome_covered_span" : "Fraction Covered",
                     "amplicons_detected" : "Amplicons Detected" }
        
        for row in reader:
            if len(header) == 0:
                for k in row.keys():
                    if k in name_map:
                        header.append(name_map[k])
                    else:
                        header.append(k)
                
            # reformat
            row['file'] = filename_to_sample(row['file'])
            row['amplicons_detected'] = row['amplicons_detected'].replace(",", ", ")
            rows.append(row.values())

        write_table("|c|c|c|c|c|p{5cm}|", header, rows)

def write_tree_section():
    plot_path = args.tree_figure.format(run_name=args.run_name)

    tree_figures = pdf_to_png(plot_path)

    print("\section{Sequence Variation}")

    # Plot each depth image for the negative control samples
    for tree_png in tree_figures:
        write_image(tree_png, 0.4)

parser = argparse.ArgumentParser()
parser.add_argument('--negative-control-depth-figure', type=str, default="plots/{run_name}_depth_by_position_negative_control.pdf")
parser.add_argument('--negative-control-table', type=str, default="qc_reports/{run_name}_negative_control_report.tsv")
parser.add_argument('--tree-figure', type=str, default="plots/{run_name}_tree_snps.pdf")
parser.add_argument('--negative-control-tsv', type=str, default="")
parser.add_argument('--run-name', type=str, default="", required=True)
args = parser.parse_args()

# Preamble
#print("# QC Report for %s" % (args.run_name))
#print("QC report generated on %s by ncov-tools" % (datetime.today().strftime('%Y-%m-%d')))

write_preamble()

print(r"\section{Summary}")
print("QC report generated for %s on %s by ncov-tools" % (escape_latex(args.run_name), datetime.today().strftime('%Y-%m-%d')))

write_tree_section()

write_negative_control_section()

write_postamble()


