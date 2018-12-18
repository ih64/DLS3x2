import pandas as pd
import numpy as np
import treecorr
import jpIO

class Pipe:

    def __init__(self):
        #set up the source and lens tables
        self.io = jpIO.io()
        self.source_table = self.io.setup_source()
        self.lens_table = self.io.setup_lens()
        self.randoms = self.io.read_randoms('randoms.csv')
        
        #tomography
        self.lens_weights = ['w1', 'w2', 'w3']
        self.source_weights = ['w1', 'w2', 'w3']

    def run(self):
        # do w theta correlations
        # cross correlations are symmetric 
        i = 0 
        for wi in self.lens_weights:
            j = 0
            for wj in self.lens_weights:
                if j > i:
                    continue
                # catch the auto correlation
                elif i == j:
                    corr = self.wtheta(self.lens_table, wi)
                    self.io.write_corrs('{}{}_mm.csv'.format(i, j), corr)
                # cross correlations
                else:
                    corr = self.wtheta(self.lens_table, wi, weight_key2=wj)
                    self.io.write_corrs('{}{}_mm.csv'.format(i, j), corr)
                j += 1
            i += 1

        # do shear shear correlations
        # cross correlations are symmetric
        i = 0 
        for wi in self.source_weights:
            j = 0
            for wj in self.source_weights:
                if j > i:
                    continue
                # catch the auto correlation
                elif i == j:
                    corr = self.shearshear(self.source_table, wi)
                    self.io.write_corrs('{}{}_gg.csv'.format(i, j), corr)
                # cross correlations
                else:
                    corr = self.shearshear(self.source_table, wi, weight_key2 = wj)
                    self.io.write_corrs('{}{}_gg.csv'.format(i, j), corr)
                j += 1
            i += 1
    
        # tangential shear correlations
        # not symmetric
        for wi, i in zip(self.lens_weights, range(len(self.lens_weights))):
            for wj, j in zip(self.source_weights, range(len(self.source_weights))):
                corr = self.gammat(self.lens_table, self.source_table, wi, wj)
                self.io.write_corrs('{}{}_gm.csv'.format(i, j), corr)
        return

    def wtheta(self, table, weight_key, weight_key2=None, **kwargs):        
        #setup correlation objects
        dd = treecorr.NNCorrelation(min_sep=1.0, max_sep=80, nbins=10, sep_units='arcmin')
        rand = treecorr.Catalog(ra=self.randoms['ra'].values, dec=self.randoms['dec'].values,
                                ra_units='radians', dec_units='radians')
        rr = treecorr.NNCorrelation(min_sep=1.0, max_sep=80, nbins=10, sep_units='arcmin')
        dr = treecorr.NNCorrelation(min_sep=1.0, max_sep=80, nbins=10, sep_units='arcmin')

        #deal with second catalog if need be
        if weight_key2 is not None:
            cat = self.io.df_to_corr(table, weight_key)
            cat2 = self.io.df_to_corr(table, weight_key2)

            rd = treecorr.NNCorrelation(min_sep=1.0, max_sep=80,
                                        nbins=10, sep_units='arcmin')

            rr.process(rand)
            dd.process(cat, cat2)
            dr.process(cat, rand)
            rd.process(rand, cat2)

            xi, varxi = dd.calculateXi(rr, dr, rd)
            sig = np.sqrt(varxi)
            r = np.exp(dd.meanlogr)
#            Coffset = calcC(rr)
            return {"xi":xi, "sig":sig, "r":r}

        #otherwise just deal with the auto correlation
        else:
            cat = self.io.df_to_corr(table, weight_key)
            #calculate w of theta given our sanitized randoms and catalog data
            rr.process(rand)
            dd.process(cat)
            dr.process(cat, rand)

            xi, varxi = dd.calculateXi(rr, dr)
            sig = np.sqrt(varxi)
            r = np.exp(dd.meanlogr)
#            Coffset = calcC(rr)
            return {"xi":xi, "sig":sig, "r":r}

    def gammat(self, lens, sources, lens_weight_key, source_weight_key, **kwargs):
        lens_corr = self.io.df_to_corr(lens, lens_weight_key)
        source_corr = self.io.df_to_corr(sources, source_weight_key, shears=True)
        rand = treecorr.Catalog(ra=self.randoms['ra'].values, dec=self.randoms['dec'].values,
                                ra_units='radians', dec_units='radians')
        #now make correlation functions
        GGL = treecorr.NGCorrelation(min_sep=0.1, max_sep=90, nbins=10, sep_units='arcmin')
        GGL.process(lens_corr, source_corr)

        # calculate random signal
        GGL_rand = treecorr.NGCorrelation(min_sep=0.1, max_sep=90, nbins=10, sep_units='arcmin')
        GGL_rand.process(rand, source_corr)
        return {'xi+': GGL.xi - GGL_rand.xi, 'xi-' : GGL.xi_im,
                'r':np.exp(GGL.meanlogr), 'sig':np.sqrt(GGL.varxi)}

    def shearshear(self, cat, weight_key1, weight_key2=None, **kwargs):
        ggkwargs = {'min_sep':1, 'max_sep':90, 'nbins':8, 'sep_units':'arcmin'}
        gg = treecorr.GGCorrelation(**ggkwargs)

        tree_cat1 = self.io.df_to_corr(cat, weight_key1, shears=True)

        if weight_key2 is not None:
            tree_cat2 = self.io.df_to_corr(cat, weight_key2, shears=True)
            gg.process(tree_cat1, tree_cat2)
        else:
            gg.process(tree_cat1)
        r = np.exp(gg.meanlogr)
        xip = gg.xip
        xim = gg.xim
        sig = np.sqrt(gg.varxi)

        return {'xip': xip, 'xim': xim, 'r':r, 'sig':sig}
    
