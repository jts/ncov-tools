'''
Suite of tests for the PerBaseCoverage module
'''

import unittest
import ncov.parser

test_coverage = ncov.parser.PerBaseCoverage(file='data/sampleA.per_base_coverage.bed')

class CoverageTest(unittest.TestCase):
    def test_get_coverage_stats(self):
        cov_stats = test_coverage.get_coverage_stats()
        expected_mean_sequencing_depth = 679.4
        expected_median_sequencing_depth = 682
        self.assertEqual(cov_stats['mean_sequencing_depth'], expected_mean_sequencing_depth)
        self.assertEqual(cov_stats['median_sequencing_depth'], expected_median_sequencing_depth)
