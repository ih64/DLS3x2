import unittest
import pandas as pd
import jpPipe


class MICETestCase(unittest.TestCase):
    ''' test suite to make sure the pipeline runs with the MICE catalog'''
    def setUp(self):
        t = pd.read_hdf('../catalogs/MICE.hdf', '/data')
        self.Pipe = jpPipe.Pipe(t)

    def tearDown(self):
        del self.Pipe

    def test_tomog_bins(self):
        ''' make sure there are tomographic bins'''
        n_lens_groups = self.Pipe.lens_groups.ngroups
        n_source_groups = self.Pipe.source_groups.ngroups
        self.assertTrue(n_lens_groups > 0, 'no lens bins')
        self.assertTrue(n_source_groups > 0, 'no source bins')
        
    def test_lens_members(self):
        ''' make sure lens bins have galaxies in them'''
        for (key, group) in self.Pipe.lens_groups:
            group_size = len(group)
            self.assertTrue(group_size > 0, 'no galaxies in lens bin')

    def test_source_members(self):
        ''' make sure the source bins have galaxies in them'''
        for (key, group) in self.Pipe.source_groups:
            group_size = len(group)
            self.assertTrue(group_size > 0, 'no galaxies in source bin')

    def test_random_density(self):
        ''' make sure the randoms are correctly scaled for each lens bin'''
        for (key, group), i in zip(self.Pipe.lens_groups, range(0, self.Pipe.lens_groups.ngroups)):
            # rand is a treecorr catalog
            rand = self.Pipe.io.read_randoms('MICE_randoms_{}.hdf'.format(i))
            group_size = len(group)
            rand_size = rand.ntot
            self.assertEqual(6*group_size, rand_size, 'randoms and lenses are not scaled correctly')

    
if __name__ == "__main__":
    unittest.main()
