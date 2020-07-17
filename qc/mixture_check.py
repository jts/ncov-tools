import argparse
import pysam
import sys
import csv
import os

from collections import defaultdict

class Variant:
    def __init__(self, contig, pos, ref, alt):
        self.contig = contig
        self.pos = pos
        self.ref = ref
        self.alt = alt
    
    def __hash__(self):
        return hash( (self.contig, self.pos, self.ref, self.alt) )
    
    def __eq__(self, other):
        return self.contig == other.contig and \
               self.pos == other.pos and \
               self.ref == other.ref and \
               self.alt == other.alt

    def __repr__(self):
        return self.to_string()

    def __str__(self):
        return self.to_string()

    def to_string(self):
        return str(self.pos) + ":" + self.ref + ">" + self.alt


class PileupObservation:
    def __init__(self, b, q):
        self.base = b
        self.quality = q

reference_contig = "MN908947.3"
#
def get_obs_at_position(samfile, position):
    observations = list()
    # -1 because pileup takes 0-based coordinates
    for pc in samfile.pileup(contig=reference_contig, start=position - 1, stop=position, truncate=True):
        for pr in pc.pileups:
            if pr.indel != 0 or pr.is_refskip or pr.query_position is None:
                continue
            else:
                b = pr.alignment.query_sequence[pr.query_position]
                q = pr.alignment.query_qualities[pr.query_position]
                po = PileupObservation(b, q)
                observations.append(po)
    return observations

def contam_frac(a, b):
    t = a + b
    c = 0
    if t > 0:
        c = b / t
    return c

def count_alleles(allele_a, allele_b, observations):
    a = 0
    b = 0

    for o in observations:
        if o.base == allele_a:
            a += 1
        elif o.base == allele_b:
            b += 1
    return (a, b)

def count_variant_support_from_bam(bam_file, variants):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    out = dict()
    for var in variants:
        observations = get_obs_at_position(samfile, var.pos)
        (count_a, count_b) = count_alleles(var.ref, var.alt, observations)
        out[var.pos] = { var.ref:count_a, var.alt:count_b }
    return out

def count_variant_support_from_fpileup(fpileup_file, variants):
    out = dict()
    variants_by_position = defaultdict(list)

    for v in variants:
        variants_by_position[v.pos].append(v)
    
    with open(fpileup_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for record in reader:
            p = int(record['position'])
            if p not in variants_by_position:
                continue

            for var in variants_by_position[p]:
                key_a = "count_" + var.ref
                key_b = "count_" + var.alt
                (count_a, count_b) = (int(record[key_a]), int(record[key_b]))
                out[var] = { var.ref:count_a, var.alt:count_b }
    return out

def calculate(sample_a, sample_b, variants, counts):
    set_a = set(variants[sample_a])
    set_b = set(variants[sample_b])

    intersection = set_a.intersection(set_b)
    symdiff = set_a.symmetric_difference(set_b)

    #debug
    #symdiff = intersection
    total_count_a = 0
    total_count_b = 0
    contaminated_variants = 0

    per_site_out = list()
    num_variants = 0
    for var in symdiff:

        # if we don't have a pileup record for this variant it didn't have
        # enough coverage, skip
        if var not in counts[sample_a]:
            continue
        
        c = counts[sample_a][var]

        # assign ref/alt to a or b sample
        allele_a = None
        allele_b = None
        if var in set_a:
            assert(var not in set_b)
            allele_a = var.alt
            allele_b = var.ref
        else:
            assert(var in set_b)
            assert(var not in set_a)
            allele_a = var.ref
            allele_b = var.alt
     

        (count_a, count_b) = (c[allele_a], c[allele_b])
        if count_a + count_b < args.min_depth:
            continue

        num_variants += 1   
        var_cf = contam_frac(count_a, count_b)
        if var_cf > args.threshold:
            contaminated_variants += 1
        per_site_out.append(str(var) + ":" + "%.3f" % (var_cf))
        total_count_a += count_a
        total_count_b += count_b

    cf = contam_frac(total_count_a, total_count_b)
    if (num_variants >= args.min_variants and contaminated_variants >= num_variants - 2) or args.show_all:
        site_str = ",".join(per_site_out)
        print("\t".join([str(x) for x in [sample_a, sample_b, num_variants, contaminated_variants, total_count_a, total_count_b, "%.3f" % (cf), site_str]]))

def load_ivar_variants(variants_fn):
    out = list()
    valid_alts = { "A", "C", "G", "T" }

    with open(variants_fn) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            alt = row['ALT']
            if alt not in valid_alts or row['PASS'] == "FALSE" or float(row['ALT_FREQ']) < args.min_allele_frequency:
                continue
            v = Variant(row['REGION'], int(row['POS']), row['REF'], row['ALT'])
            out.append(v)
    return out

def load_msa_alleles(alleles_fn):
    out = dict()
    valid_alts = { "A", "C", "G", "T" }
    with open(alleles_fn) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            sample = row['name'].replace("Consensus_", "")
            if sample not in out:
                out[sample] = list()
            ra = row['ref_allele']
            aa = row['alt_allele']
            v = Variant("contig", int(row['pos']), row['ref_allele'], row['alt_allele'])
            if aa in valid_alts:
                out[sample].append(v)
    return out
#
def load_fofn(fofn_fn):
    out = dict()
    with open(fofn_fn) as fh:
        for line in fh:
            path = line.rstrip()
            fn = os.path.basename(path)
            fields = fn.split(".")
            out[fields[0]] = path
    return out

parser = argparse.ArgumentParser()
parser.add_argument('--fpileup-fofn', type=str, default="")
parser.add_argument('--alleles-tsv', type=str, default="")
parser.add_argument('--show-variants', action='store_true')
parser.add_argument('--show-all', action='store_true')
parser.add_argument('--threshold', type=float, default=0.1)
parser.add_argument('--min-variants', type=int, default=4)
parser.add_argument('--min-depth', type=int, default=10)
parser.add_argument('--min-allele-frequency', type=float, default=0.0)
parser.add_argument('--sample', type=str, default="")

args = parser.parse_args()

# read fofn files
# each file must be named such that the sample is the first part of
# the filename
fpileup_files = load_fofn(args.fpileup_fofn)

# load variants observed in the MSA
variants = load_msa_alleles(args.alleles_tsv)

# make a set of all variants
variant_union = set()
for k,v in variants.items():
    variant_union = variant_union.union(v)

# count bases at each variant position for each sample
all_samples = variants.keys()
counts = dict()
for s in all_samples:

    if s in fpileup_files:
        counts[s] = count_variant_support_from_fpileup(fpileup_files[s], variant_union)
    else:
        sys.stderr.write("Could not find sequencing data for %s\n" % (s))
        sys.exit(1)

# determine sample set to use
outer_samples = all_samples
if args.sample != "":
    outer_samples = [args.sample]

print("\t".join(["sample_a", "sample_b", "variants_checked", "variants_mixed", \
                 "read_support_allele_a", "read_support_allele_b", "mixture_fraction", \
                 "mixed_sites"]))

for a_sample in outer_samples:
    for b_sample in all_samples:
        if a_sample == b_sample:
            continue
        calculate(a_sample, b_sample, variants, counts)
