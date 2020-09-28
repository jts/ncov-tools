import csv
import sys
import os

print("\t".join(["sample", "num_subs", "num_c_to_t",  "num_indels"]))

for fn in sys.argv[1:]:
    with open(fn) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        sample_name = os.path.basename(fn).split(".")[0]
        num_subs = 0
        num_ct = 0
        num_indels = 0
        for row in reader:
            if row['ALT'][0] == "+" or row['ALT'][0] == '-':
                num_indels += 1
            else:
                num_subs += 1
                if row['REF'] == 'C' and row['ALT'] == 'T':
                    num_ct += 1
        print("\t".join([str(x) for x in [sample_name, num_subs, num_ct, num_indels]]))
