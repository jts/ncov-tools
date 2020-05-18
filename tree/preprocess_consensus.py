import collections
import pysam
import sys

completeness_threshold = 0.75

for fn in sys.argv[1:]:
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
