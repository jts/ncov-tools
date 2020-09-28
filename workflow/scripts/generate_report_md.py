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

def write_image(image_filename):
    print("![img](%s)" % (image_filename))

def write_table(header, width_map, rows):
    print("")
    print("|".join(header))

    span = list()
    for h in header:
        w = 5
        if h in width_map:
            w = width_map[h]
        span.append("-" * w)
        
    print("|".join(span))
    for row in rows:
        print("|".join(row))
    print("")

def filename_to_sample(fn):
    return os.path.basename(fn).split(".")[0]

def write_negative_control_section():

    plot_path = args.negative_control_depth_figure.format(run_name=args.run_name)

    nc_depth_figures = pdf_to_png(plot_path)

    print("\n## Negative Control Report")

    # Plot each depth image for the negative control samples
    for depth_png in nc_depth_figures:
        write_image(depth_png)

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
        
        width_map = { "Covered" : 2,
                      "qc": 2 }
        
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

        write_table(header, width_map, rows)

def write_tree_section():
    plot_path = args.tree_figure.format(run_name=args.run_name)

    tree_figures = pdf_to_png(plot_path)

    print("\n## Sequence Variation Tree")

    # Plot each depth image for the negative control samples
    for tree_png in tree_figures:
        write_image(tree_png)

parser = argparse.ArgumentParser()
parser.add_argument('--negative-control-depth-figure', type=str, default="plots/{run_name}_depth_by_position_negative_control.pdf")
parser.add_argument('--negative-control-table', type=str, default="qc_reports/{run_name}_negative_control_report.tsv")
parser.add_argument('--tree-figure', type=str, default="plots/{run_name}_tree_snps.pdf")
parser.add_argument('--negative-control-tsv', type=str, default="")
parser.add_argument('--run-name', type=str, default="", required=True)
args = parser.parse_args()

# Preamble
print("# QC Report for %s" % (args.run_name))
print("QC report generated on %s by ncov-tools" % (datetime.today().strftime('%Y-%m-%d')))

write_tree_section()

write_negative_control_section()
