import argparse
import sys
import csv
import os

parser = argparse.ArgumentParser()
parser.add_argument('--primer-snp-bed', type=str, required=True)
parser.add_argument('--amplicon-depth-tsv', type=str, required=True)
parser.add_argument('--sample-name', type=str, default="none")
args, files = parser.parse_known_args()

def main():
    print("\t".join(["sample", "contig", "position", "ref", "alt", "primer_name", "amplicon", "depth"]))
    
    # read in the depth of each amplicon
    depth_map = dict()
    with open(args.amplicon_depth_tsv, 'r') as ifh:
        reader = csv.DictReader(ifh, delimiter='\t')
        for record in reader:
            depth_map[record['amplicon_id']] = record['mean_depth']
    
    # read in the bed file containing mutations that overlap primers
    seen = dict()
    with(open(args.primer_snp_bed, 'r')) as ifh:
        for line in ifh:
            fields = line.rstrip().split()
            assert(len(fields) == 16)
            primer_name = fields[13]
            primer_direction_idx = max(primer_name.find("_LEFT"), primer_name.find("_RIGHT"))
            if primer_direction_idx == -1:
                sys.stderr.write("Could not parse primer name %s" % (primer_name))
                sys.exit(1)
            amplicon_id = primer_name[0:primer_direction_idx]

            # skip duplicate entries when a mutation lands in ALT versions of primers
            key = ":".join([fields[0], fields[1], fields[3], fields[4], amplicon_id])
            if key not in seen:
                print("\t".join([args.sample_name, fields[0], fields[1], fields[3], fields[4], primer_name, amplicon_id, depth_map[amplicon_id]]))
                seen[key] = 1

if __name__ == '__main__':
    main()
