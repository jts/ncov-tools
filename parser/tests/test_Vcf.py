'''
Suite of tests for the Vcf module
'''

import unittest
import ncov.parser

test_variants = ncov.parser.Vcf(file='data/sampleA.pass.vcf')

class VariantsTest(unittest.TestCase):
    def test_get_variant_counts(self):
        var_counts = test_variants.get_variant_counts()
        self.assertEqual(var_counts['num_variants_snvs'], 5)
        self.assertEqual(var_counts['num_variants_indel'], 3)
        self.assertEqual(var_counts['num_variants_indel_triplet'], 1)
