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

# process an fpileup.tsv file into a dictionary
# of observed read counts for each variant in the union
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

# if the sample has an entry for the input var, return the consensus allele
# otherwise return reference
def assign_allele(var, sample, variants):
    if var in variants[sample]:
        return variants[sample][var]
    else:
        return var.ref

# check a pair of samples for evidence of contamination/mixtures
def calculate(sample_a, sample_b, variants, counts):

    union = dict()
    for v in variants[sample_a]:
        union[v] = 1
    for v in variants[sample_b]:
        union[v] = 1

    # output stats
    total_count_a = 0
    total_count_b = 0
    contaminated_variants = 0
    per_site_out = list()
    num_variants = 0
    alleles_by_position = dict()
    for var in union:
        # assign alleles to each sample
        allele_b = assign_allele(var, sample_b, variants)
        allele_a = assign_allele(var, sample_a, variants)

        # if B has an IUPAC code it is impossible to resolve - skip this position
        if allele_b not in "ACGT":
            continue

        # if A is an IUPAC code, resolve to the non-B allele
        allele_a = resolve_iupac(allele_a, allele_b)

        # skip non-informative positions
        if allele_a == allele_b:
            continue

        # skip the sample if it isn't in counts
        if sample_a not in counts:
            continue

        # skip positions without any depth
        if var not in counts[sample_a]:
            continue

        c = counts[sample_a][var]
        # look up counts for a and b allele, check depth
        (count_a, count_b) = (c[allele_a], c[allele_b])

        if count_a + count_b < args.min_depth:
            continue

        # informative variant, count it
        num_variants += 1

        # calculate mixture proportion
        var_cf = contam_frac(count_a, count_b)
        if var_cf > args.threshold:
            contaminated_variants += 1
        #print(var.pos, allele_a, count_a, allele_b, count_b, var_cf)

        per_site_out.append(str(var) + ":" + "%.3f" % (var_cf))
        total_count_a += count_a
        total_count_b += count_b

    # output if sufficient evidence of mixture
    cf = contam_frac(total_count_a, total_count_b)
    if (num_variants >= args.min_variants and contaminated_variants >= num_variants - 2) or args.show_all:
        site_str = ",".join(per_site_out)
        print("\t".join([str(x) for x in [sample_a, sample_b, num_variants, contaminated_variants, total_count_a, total_count_b, "%.3f" % (cf), site_str]]))

# load variants from the ivar variants.tsv file
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

iupac_map = { "M": ["A", "C"],
              "R": ["A", "G"],
              "W": ["A", "T"],
              "S": ["C", "G"],
              "Y": ["C", "T"],
              "K": ["G", "T"] }

# given an IUPAC code (like Y) and a fixed other base (like C) return the other base (T)
def resolve_iupac(code, fixed_base):
    if code in ["A", "C", "G", "T"]:
        return code # not ambiguous
    if code in iupac_map:
        (a, b) = iupac_map[code]
        if a == fixed_base:
            return b
        elif b == fixed_base:
            return a
        else:
            assert(False)
            return None
    else:
        return None
        # triallelic
        assert(False)

# load observed variants from the output of align2alleles.py
def load_msa_alleles(alleles_fn):
    out = dict()
    with open(alleles_fn) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            sample = row['name'].replace("Consensus_", "")
            if sample not in out:
                out[sample] = dict()
            ra = row['ref_allele']
            ca = row['alt_allele'] # consensus allele, can be ambiguous
            
            # ignore Ns and singletons
            if ca == "N" or int(row['samples_with_allele']) == 1:
                continue
            aa = resolve_iupac(ca, ra)

            v = Variant("contig", int(row['pos']), ra, aa)
            out[sample][v] = ca
    return out

# load a map from sample name -> file path from a file-of-filenames (fofn)
def load_fofn(fofn_fn):
    out = dict()
    with open(fofn_fn) as fh:
        for line in fh:
            path = line.rstrip()
            fn = os.path.basename(path)
            fields = fn.split(".")
            out[fields[0]] = path
    return out

# 
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

# count variants and remove singletons
var_count = defaultdict(int)
for sample in variants:
    for var in variants[sample]:
        var_count[var] += 1

# create the union of variants
variant_union = set()
for var, count in var_count.items():
    if count >= 1:
        variant_union.add(var)

# count bases at each variant position for each sample
all_samples = variants.keys()
counts = dict()
for s in all_samples:

    if s in fpileup_files:
        counts[s] = count_variant_support_from_fpileup(fpileup_files[s], variant_union)
    else:
        sys.stderr.write("Could not find sequencing data for %s\n" % (s))
        #sys.exit(1)

# determine sample set to use
outer_samples = all_samples
if args.sample != "":
    outer_samples = [args.sample]

# output header
print("\t".join(["sample_a", "sample_b", "variants_checked", "variants_mixed", \
                 "read_support_allele_a", "read_support_allele_b", "mixture_fraction", \
                 "mixed_sites"]))

# check all pairs of samples for contamination
for a_sample in outer_samples:
    for b_sample in all_samples:
        if a_sample == b_sample:
            continue
        calculate(a_sample, b_sample, variants, counts)
