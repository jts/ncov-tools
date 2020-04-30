import pysam
import sys

for fn in sys.argv[1:]:
    fa = pysam.FastxFile(fn)
    for record in fa:

        # shrink the long name ivar uses
        record.name = record.name.replace(".primertrimmed.consensus_threshold_0.75_quality_20", "")
        sys.stdout.write(">" + record.name + "\n")
        sys.stdout.write(record.sequence + "\n")
