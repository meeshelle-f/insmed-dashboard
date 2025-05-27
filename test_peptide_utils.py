import unittest 
from streamlit_app_v2 import calculate_peptide_similarity, determine_most_similar_groups

class TestCalculatePeptideSimilarity(unittest.TestCase):

    def test_exact_match(self):
        self.assertEqual(calculate_peptide_similarity("ABC", "ABC"), 1.0)

    def test_no_match(self):
        self.assertEqual(calculate_peptide_similarity("ABC", "XYZ"),  0.0)

    def test_partial_match(self):
        self.assertAlmostEqual(calculate_peptide_similarity("ABCDEF", "ABXCDE"), 2/6)
        self.assertAlmostEqual(calculate_peptide_similarity("HELLO", "HELBO"), 4/5)

    def test_empty_peptide(self):
        self.assertEqual(calculate_peptide_similarity("", "ABC"), 0.0)
        self.assertEqual(calculate_peptide_similarity("ABC", ""), 0.0)
        self.assertEqual(calculate_peptide_similarity("", ""), 0.0)

    def test_whitespace_handling(self):
        self.assertEqual(calculate_peptide_similarity(" ABC ", "ABC"), 1.0) 
        self.assertEqual(calculate_peptide_similarity("\nABC\t", "ABC"), 1.0)

    def test_different_lengths(self):
        self.assertEqual(calculate_peptide_similarity("ABC", "AB"), 0.0)

    def test_non_string_input(self):
        self.assertEqual(calculate_peptide_similarity(123, "123"), 1.0)
        self.assertEqual(calculate_peptide_similarity(None, "ABC"), 0.0)
        self.assertEqual(calculate_peptide_similarity("ABC", 123), 0.0)

if __name__ == '__main__':
    unittest.main()

class TestDetermineMostSimilarGroups(unittest.TestCase):

    def test_normal_case(self):
        self.df = pd.DataFrame({'group': ['G1', 'G1', 'G2', 'G2'],'peptide': ['AB', 'AB', 'AA', 'BB']})
        result = determine_most_similar_groups(self.df)
        self.assertEqual(result, {'group_pair': ('G1', 'G2'), 'average_similarity_score': 0.5, 'all_pair_scores': {('G1', 'G2'): 0.5}})

if __name__ == '__main__':
    unittest.main()
