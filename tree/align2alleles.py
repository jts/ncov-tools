import pysam
import sys
import argparse
import collections

# write results in matrix form where rows are samples
# and columns are variant positions
def write_result_matrix(variant_positions, sequences):
    # print header
    # +1 is to report 1-based coordinates
    print("\t".join(["strain"] + [str(i + 1) for i in variant_positions]))

    for s in sequences:
        out = list()
        out.append(s[0])

        for i in variant_positions:
            out.append(str(s[1][i]))

        print("\t".join(out))

# write results as a TSV file with one row per variant found in a sample
def write_variant_list(variant_positions, sequences, reference_sample):
    print("\t".join(["name", "pos", "ref_allele", "alt_allele"]))

    for s in sequences:
        for i in variant_positions:
            b = s.sequence[i]
            if b != reference_sample.sequence[i]:
                print("\t".join([s.name, str(i+1), reference_sample.sequence[i], b]))

#
#
#
parser = argparse.ArgumentParser()
parser.add_argument('--min-allele-count', default=1, type=int)
parser.add_argument('--mode', default="variant_list", type=str)
parser.add_argument('--reference-name', default="Wuhan-Hu-1/2019", type=str)
args, extra = parser.parse_known_args()

fa = pysam.FastxFile(extra[0])
sequences = list()

# read all sequence from the MSA
alen = -1
for r in fa:
    if alen == -1:
        alen = len(r.sequence)
    r.sequence = r.sequence.upper()
    sequences.append(r)

# discover variant columns
counters = list()
for i in range(0, alen):
    counters.append(collections.Counter())

#
reference_sample = None

# count the number of occurrences of each base at each position
for s in sequences:

    if s.name == args.reference_name:
        reference_sample = s

    for i in range(0, alen):
        b = s.sequence[i]
        counters[i].update(b)

if reference_sample is None:
    sys.stderr.write("error: reference sample could not be found")
    sys.exit(1)

# write results
if args.mode == "variant_list":
    # determine which positions have a SNP
    variant_positions = list()
    for i in range(0, alen):
        supported_bases = list()
        for b in ("A", "C", "G", "T"):
            c = counters[i][b]
            if c >= args.min_allele_count:
                supported_bases.append(b)

        if len(supported_bases) > 1:
            variant_positions.append(i)
    write_variant_list(variant_positions, sequences, reference_sample)

elif args.mode == "variant_frequency":
    for i in range(0, alen):
        reference_base = reference_sample.sequence[i]
        for b in ("A", "C", "G", "T"):
            if b == reference_base:
                continue
            c = counters[i][b]
            if c > 0:
                out = [ reference_name, str(i + 1), reference_base, b, counters[i][reference_base], counters[i][b] ]
                print("\t".join([str(x) for x in out]))
