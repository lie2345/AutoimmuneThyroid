import importlib.util
import math
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('stats',Path(__file__).with_name('57_run_external_validation.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
C={'planned_tests':20,'family_alpha':.05,'equivalence_or_upper':1.10,'precision_power_threshold':.8}


class Statistics(unittest.TestCase):
    def test_direction_and_source_p(self):
        a=m.calculate(.2,.05,.001,.3,.03,C)
        self.assertEqual(a['classification'],'same_direction_association_support')
        b=m.calculate(-.2,.05,.001,.3,.03,C)
        self.assertEqual(b['classification'],'opposite_direction_association_support')
        self.assertNotEqual(a['normal_p'],.001)
        self.assertEqual(m.calculate(.2,.05,.1,.3,.03,C)['classification'],'inconclusive_low_precision_OR1_10_scenario')

    def test_equivalence_independent_of_association(self):
        a=m.calculate(.01,.001,1e-10,.2,.05,C)
        self.assertTrue(a['equivalence_flag'])
        self.assertEqual(a['classification'],'same_direction_association_support')
        self.assertGreater(a['equivalence_beta_upper995']-a['equivalence_beta_lower995'],2*1.96*.001)
        self.assertEqual(m.calculate(0,.01,1,.2,.05,C)['classification'],'no_primary_support_small_effect_equivalence')

    def test_heterogeneity_and_symmetry(self):
        a=m.calculate(.2,.1,.05,.2,.1,C)
        self.assertEqual(a['heterogeneity_q'],0)
        self.assertEqual(a['heterogeneity_p'],1)
        self.assertEqual(a['i2_descriptive'],0)
        b=m.calculate(.2,.1,.05,-.2,.1,C)
        self.assertAlmostEqual(b['heterogeneity_q'],8)
        self.assertAlmostEqual(b['i2_descriptive'],.875)
        c=m.calculate(-.2,.1,.05,.2,.1,C)
        self.assertEqual(c['heterogeneity_p'],b['heterogeneity_p'])

    def test_power_and_mde(self):
        self.assertAlmostEqual(m.precision_power(0,.1,.0025),.0025)
        self.assertAlmostEqual(m.precision_power(math.log(m.mde(.1,.0025)),.1,.0025),.8)
        self.assertGreater(m.precision_power(.2,.1,.0025),m.precision_power(.1,.1,.0025))
        with self.assertRaises(AssertionError):m.calculate(0,0,.1,.1,.1,C)


if __name__=='__main__':unittest.main()
