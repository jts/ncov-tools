import argparse
import numpy
import sys
import csv
import os

from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('--input-tsv', type=str, required=True)
args, files = parser.parse_known_args()

def main():
    print("\t".join(["contig", "position", "ref", "alt", "primer_name", "amplicon", "num_samples", "mean_depth", "min_depth", "max_depth"]))
    
    # read in the depth per sample
    depth_map = defaultdict(list)

    with open(args.input_tsv, 'r') as ifh:
        reader = csv.DictReader(ifh, delimiter='\t')
        for record in reader:
            key = ":".join([record['contig'], record['position'], record['ref'], record['alt'], record['primer_name'], record['amplicon']])
            depth_map[key].append(float(record['depth']))
    

    for key in depth_map.keys():
        mean_depth = numpy.mean(depth_map[key])
        min_depth = numpy.min(depth_map[key])
        max_depth = numpy.max(depth_map[key])
        fields = key.split(":")
        fields.append(str(len(depth_map[key])))
        fields.append("%.2f" % mean_depth)
        fields.append("%.2f" % min_depth)
        fields.append("%.2f" % max_depth)
        print("\t".join(fields))

if __name__ == '__main__':
    main()
