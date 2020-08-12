'''
Suite of tests for the Variants module
'''

import unittest
import ncov.parser

test_variants = ncov.parser.Variants(file='data/sampleA.variants.tsv')

class VaiantsTest(unittest.TestCase):
    def test_get_total_variants(self):
        var_counts = test_variants.get_total_variants()
        self.assertEqual(var_counts['num_variants_snvs'], 9)
        self.assertEqual(var_counts['num_variants_indel'], 1)
        self.assertEqual(var_counts['num_variants_indel_triplet'], 1)
