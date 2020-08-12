'''
Suite of tests for the Consensus module
'''

import unittest
import ncov.parser

test_consensus = ncov.parser.Consensus(file='data/tester.fa')


class ConsensusTest(unittest.TestCase):
    def test_count_iupac_in_fasta(self):
        fasta_iupac_counts = test_consensus.count_iupac_in_fasta()
        self.assertEqual(fasta_iupac_counts['num_consensus_iupac'], 9)
        self.assertEqual(fasta_iupac_counts['num_consensus_n'], 13)

    def test_get_genome_completeness(self):
        fasta_iupac_counts = test_consensus.get_genome_completeness()
        self.assertEqual(fasta_iupac_counts['genome_completeness'], 0.8992)