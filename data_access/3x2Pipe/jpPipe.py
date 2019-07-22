import os
import argparse
from os.path import join
import pandas as pd
import numpy as np
import treecorr
import jpIO

class Pipe:

    def __init__(self, infile, filename='config.yaml', 
                 inpath='/global/cscratch1/sd/ihasan/catalogs', outdir='/global/cscratch1/sd/ihasan/flask/output'):
        #set up the source and lens tables
        self.io = jpIO.io(filename)
        #TODO
        #use flask versionof settup
        fname = infile.split('.dat')[0]
        self.outdir = join(outdir,fname)
        if os.path.exists(self.outdir):
          pass
        else:
            os.mkdir(self.outdir)

        tables = self.io.setup_FLASK_cats(join(inpath, infile))
        self.source_table = tables['source']
        self.lens_table = tables['lens']
        self.randoms = tables['rand']
        #tomography
        self.lens_groups = self.lens_table.groupby('z_bin')
        self.source_groups = self.source_table.groupby('z_bin')

    def run(self):
        # do w theta correlations
        # cross correlations are symmetric
        for (keyi, groupi), i in zip(self.lens_groups, range(0, 3)):
            for (keyj, groupj), j in zip(self.lens_groups, range(0, 3)):
                if (j > i):
                    continue
                # catch the auto correlation
                elif i == j:
                    corr = self.wtheta(groupi, i)
                    self.io.write_corrs('{}{}_mm.csv'.format(i, j), corr, self.outdir)
                # cross correlations
                else:
                    corr = self.wtheta(groupi, i, table2=groupj, bin_number_2=j)
                    self.io.write_corrs('{}{}_mm.csv'.format(i, j), corr, self.outdir)

        # do shear shear correlations
        # cross correlations are symmetric
        for (keyi, groupi), i in zip(self.source_groups, range(0, 3)):
            for (keyj, groupj), j in zip(self.source_groups, range(0, 3)):
                if (j > i):
                    continue
                # catch the auto correlation
                elif i == j:
                    corr = self.shearshear(groupi)
                    self.io.write_corrs('{}{}_gg.csv'.format(i, j), corr, self.outdir)
                # cross correlations
                else:
                    corr = self.shearshear(groupi, cat2=groupj)
                    self.io.write_corrs('{}{}_gg.csv'.format(i, j), corr, self.outdir)
    
        # tangential shear correlations
        # not symmetric
        for (keyl, groupl), i in zip(self.lens_groups, range(0, 3)):
            for (keys, groups), j in zip(self.source_groups, range(0, 3)):
                corr = self.gammat(groupl, groups, i)
                self.io.write_corrs('{}{}_gm.csv'.format(i, j), corr, self.outdir)
        return

    def wtheta(self, table, bin_number, table2=None, bin_number_2=None):
        '''calculate position position correlation'''
        #setup correlation objects, random catalog
        corr_kwargs = {'min_sep':1.0, 'max_sep':90, 'nbins':10,
                       'sep_units':'arcmin'}
        dd = treecorr.NNCorrelation(**corr_kwargs)
        rr = treecorr.NNCorrelation(**corr_kwargs)
        dr = treecorr.NNCorrelation(**corr_kwargs)

        rand = self.randoms

#        assert len(table)*6 == rand.ntot, "randoms are not scaled correctly for auto"
        
        #deal with second catalog if need be
        if table2 is not None:
            cat = self.io.df_to_corr(table)
            cat2 = self.io.df_to_corr(table2)
            rand2 = self.randoms

#            assert len(table2)*6 == rand2.ntot, "randoms are not scaled correctly for cross"
            
            rd = treecorr.NNCorrelation(**corr_kwargs)

            rr.process(rand, rand2)
            dd.process(cat, cat2)
            dr.process(cat, rand2)
            rd.process(rand, cat2)

            xi, varxi = dd.calculateXi(rr, dr, rd)
            sig = np.sqrt(varxi)
            r = np.exp(dd.meanlogr)
#            Coffset = calcC(rr)
            return {"xi":xi, "sig":sig, "r":r}

        #otherwise just deal with the auto correlation
        else:
            cat = self.io.df_to_corr(table)
            #calculate w of theta given our sanitized randoms and catalog data
            rr.process(rand)
            dd.process(cat)
            dr.process(cat, rand)

            xi, varxi = dd.calculateXi(rr, dr)
            sig = np.sqrt(varxi)
            r = np.exp(dd.meanlogr)
#            Coffset = calcC(rr)
            return {"xi":xi, "sig":sig, "r":r}

    def gammat(self, lens, sources, lens_bin_idx):
        '''calculate tangential shear correlation'''
        lens_corr = self.io.df_to_corr(lens, shears=False)
        source_corr = self.io.df_to_corr(sources, shears=True)
        rand = self.randoms

        corr_kwargs = {'min_sep':3, 'max_sep':90, 'nbins':6, 'sep_units':'arcmin'}
        #now make correlation functions
        GGL = treecorr.NGCorrelation(**corr_kwargs)
        GGL.process(lens_corr, source_corr)

        # calculate random signal
        GGL_rand = treecorr.NGCorrelation(**corr_kwargs)
        GGL_rand.process(rand, source_corr)
        return {'xi+': GGL.xi - GGL_rand.xi, 'xi-' : GGL.xi_im,
                'r':np.exp(GGL.meanlogr), 'sig':np.sqrt(GGL.varxi)}

    def shearshear(self, cat1, cat2=None):
        '''calculate shear-shear correlation '''
        ggkwargs = {'min_sep':1, 'max_sep':90, 'nbins':8, 'sep_units':'arcmin'}
        gg = treecorr.GGCorrelation(**ggkwargs)
        tree_cat1 = self.io.df_to_corr(cat1, shears=True)
        if cat2 is not None:
            tree_cat2 = self.io.df_to_corr(cat2, shears=True)
            gg.process(tree_cat1, tree_cat2)
        else:
            gg.process(tree_cat1)
        r = np.exp(gg.meanlogr)
        xip = gg.xip
        xim = gg.xim
        sig = np.sqrt(gg.varxip)
        return {'xip': xip, 'xim': xim, 'r':r, 'sig':sig}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='the path to the input flask catalog')
    parser.add_argument('--outfile', help='the output path for the correlations')
    args = parser.parse_args()
    if args.outfile:
        p = Pipe(args.infile, outdir=args.outfile)
    else:
        p = Pipe(args.infile)
    p.run()
