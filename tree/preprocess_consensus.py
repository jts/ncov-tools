import collections
import pysam
import sys
import argparse

parser = argparse.ArgumentParser(description='a consensus sequence pre-processing script')
parser.add_argument('-c', '--completeness', default=0.75,
                    help='the minimum completeness threshold (default: 0.75)')
parser.add_argument('files', action='store', nargs='*')
if len(sys.argv) < 2:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

completeness_threshold = float(args.completeness)

for fn in args.files:
    fa = pysam.FastxFile(fn)
    for record in fa:

        # shrink the long name ivar uses
        record.name = record.name.replace(".primertrimmed.consensus_threshold_0.75_quality_20", "")

        # require >75% completeness (non-N)
        base_counter = collections.Counter()
        for b in record.sequence:
            base_counter.update(b.upper())

        total = 0
        for base, count in base_counter.items():
            total += count
        completeness = 0
        if total > 0:
            completeness = 1 - (float(base_counter['N']) / total)
        if completeness >= completeness_threshold:
            sys.stdout.write(">" + record.name + "\n")
            sys.stdout.write(record.sequence + "\n")
