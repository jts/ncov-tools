import argparse
import pysam
import sys
import csv
import os

class RegionCoverage:
    def __init__(self, name, threshold):
        self.name = name
        self.total_positions = 0
        self.covered_positions = 0
        self.threshold = threshold

    def update(self, position, depth):
        self.total_positions += 1
        if depth >= self.threshold:
            self.covered_positions += 1

    def span(self):
        return float(self.covered_positions) / float(self.total_positions)

    def print_tsv(self):
        s = "\t".join([str(self.name), str(self.covered_positions), str(self.total_positions), "%.3f" % (self.span())])
        print(s)

parser = argparse.ArgumentParser()
parser.add_argument('--max-coverage', type=int, default=5)
parser.add_argument('--max-amplicon-span', type=int, default=0.5)
args, files = parser.parse_known_args()


print("\t".join(["file", "qc", "genome_covered_bases", "genome_total_bases", "genome_covered_span", "amplicons_detected"]))

for fn in files:
    with open(fn) as f:
        genome_coverage = RegionCoverage("whole_genome", args.max_coverage)
        amplicon_coverage = dict()

        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            position = int(row["start"]) + int(row["position"])
            depth = int(row["depth"])
            amplicon_id = int(row["amplicon_id"])

            #
            genome_coverage.update(position, depth)
            if amplicon_id not in amplicon_coverage:
                amplicon_coverage[amplicon_id] = RegionCoverage(amplicon_id, args.max_coverage)
            amplicon_coverage[amplicon_id].update(position, depth)
         
        detected_amplicons = list()       
        for amplicon_id in sorted(amplicon_coverage.keys()):
            ac = amplicon_coverage[amplicon_id]
            if ac.span() > args.max_amplicon_span:
                detected_amplicons.append(str(ac.name))

        qc = "PASS"
        if len(detected_amplicons) > 0:
            qc = "WARN"
        if genome_coverage.span() > 0.01:
            qc = "FAIL"
        print("\t".join([fn, qc, str(genome_coverage.covered_positions), str(genome_coverage.total_positions), \
                         "%.3f" % (genome_coverage.span()), ",".join(detected_amplicons)]))
