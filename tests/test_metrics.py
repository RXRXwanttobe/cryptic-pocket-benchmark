"""Small scientific edge cases, rather than implementation snapshots."""
import unittest
from pocketbench.metrics import iou,top_hit,score_target,wilson,paired_comparison

class MetricsTests(unittest.TestCase):
    def test_top_k_is_not_union_of_pockets(self):
        # Each pocket is individually below threshold; their union would pass.
        pockets=[{"rank":i+1,"residues":[str(i)]} for i in range(5)]
        self.assertFalse(top_hit(pockets,{str(i) for i in range(5)},5,.25))

    def test_exact_threshold_and_rank_boundary(self):
        pockets=[{"rank":6,"residues":["1"]}]
        reference={"1","2","3","4"}
        self.assertFalse(top_hit(pockets,reference,5))
        self.assertTrue(top_hit(pockets,reference,6))
        self.assertEqual(iou(["1"],reference),.25)

    def test_extra_residues_penalize_iou(self):
        self.assertLess(iou(["1","2","3","4"],["1"]),iou(["1"],["1"]))

    def test_unmapped_and_duplicate_ranks_rejected(self):
        with self.assertRaises(ValueError):score_target([{"rank":1,"residues":["99"]}],{"1"},{"1"})
        with self.assertRaises(ValueError):score_target([{"rank":1,"residues":["1"]}]*2,{"1"},{"1"})

    def test_zero_hits_have_nonzero_wilson_upper_limit(self):
        lower,upper=wilson(0,100)
        self.assertAlmostEqual(lower,0)
        self.assertGreater(upper,0)
        self.assertLess(upper,.04)

    def test_paired_test_uses_discordant_targets(self):
        a=[1]*7+[0]*93;b=[0]*100
        result=paired_comparison(a,b)
        self.assertEqual(result["a_only"],7)
        self.assertEqual(result["b_only"],0)
        self.assertAlmostEqual(result["mcnemar_p"],.015625)

if __name__=="__main__":unittest.main()
