'''
Suite of tests for the Alelles module
'''

import unittest
import ncov.parser

test_alleles = ncov.parser.Alleles(file='data/alleles.tsv')


class AllelesTest(unittest.TestCase):
    def test_init(self):
        self.assertEqual(test_alleles.data['sampleA']['241']['alt'], 'T')
    
    def test_get_variant_counts(self):
        self.assertEqual(test_alleles.get_variant_counts(sample='sampleA')['num_consensus_snvs'], 5)
        self.assertEqual(test_alleles.get_variant_counts(sample='sampleB')['num_consensus_snvs'], 6)
