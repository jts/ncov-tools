import argparse
import pysam
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--reference-genome', type=str, default="", required=True)
parser.add_argument('--bam', type=str, default="", required=True)

args = parser.parse_args()
samfile = pysam.AlignmentFile(args.bam, "rb")
reference = pysam.FastaFile(args.reference_genome)

print("\t".join(["contig", "position", "depth", "ref_base", "count_A", "count_C", "count_G", "count_T", "count_del", "alt_frequency"]))
for pc in samfile.pileup(ignore_orphans=False):
    freqs = {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0, 'R': 0 }
    for pr in pc.pileups:
        if pr.is_del:
            freqs['-'] += 1
        elif pr.is_refskip:
            freqs['R'] += 1
        else:
            b = pr.alignment.query_sequence[pr.query_position]
            freqs[b] += 1

    reference_base = reference.fetch(pc.reference_name, pc.reference_pos, pc.reference_pos+1)

    # count alt depth
    bc_depth = 0
    alt_depth = 0
    for b in "ACGT":
        bc_depth += freqs[b]
        if b != reference_base:
            alt_depth += freqs[b] 

    if bc_depth > 0:
        alt_freq = alt_depth / bc_depth
    else:
        alt_freq = 0.0

    out = list()
    out.append(pc.reference_name)
    # output 1-based coordinates for consistency with samtools mpileup
    out.append(pc.pos + 1)
    out.append(bc_depth + freqs['-'])
    out.append(reference_base)
    out.append(freqs['A'])
    out.append(freqs['C'])
    out.append(freqs['G'])
    out.append(freqs['T'])
    out.append(freqs['-'])
    out.append(alt_freq)
    print("\t".join([str(x) for x in out]))
