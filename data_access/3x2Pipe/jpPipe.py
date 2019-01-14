import pandas as pd
import numpy as np
import treecorr
import jpIO

class Pipe:

    def __init__(self, full_table):
        #set up the source and lens tables
        self.io = jpIO.io()
        self.source_table = self.io.setup_source(full_table)
        self.lens_table = self.io.setup_lens(full_table)
        self.randoms = self.io.read_randoms('buzzard_randoms_sq.hdf')

        #tomography
        s_bins = pd.cut(self.source_table['z_b'], [.4, .6, .8, 1.])
        l_bins = pd.cut(self.lens_table['z_b'], [.3, .45, .6, .8])

        self.lens_groups = self.lens_table.groupby(l_bins)
        self.source_groups = self.source_table.groupby(s_bins)

    def run(self):
        # do w theta correlations
        # cross correlations are symmetric
        for (keyi, groupi), i in zip(self.lens_groups, range(0, 3)):
            for (keyj, groupj), j in zip(self.lens_groups, range(0, 3)):
                if (j > i):
                    continue
                # catch the auto correlation
                elif i == j:
                    corr = self.pos_pos_corr((groupi['alpha'], groupi['delta']), 
                                             (groupi['alpha'], groupi['delta']), 
                                             (self.randoms['ra'], self.randoms['dec']), 
                                             (self.randoms['ra'], self.randoms['dec']), 
                                             same_cell=True, same_zshell=True)
                    self.io.write_corrs('{}{}_mm.csv'.format(i, j), corr)
                # cross correlations
                else:
                    corr = self.pos_pos_corr((groupi['alpha'], groupi['delta']),
                                             (groupj['alpha'], groupj['delta']),
                                             (self.randoms['ra'], self.randoms['dec']),
                                             (self.randoms['ra'], self.randoms['dec']), 
                                             same_cell=True, same_zshell=False)
                    self.io.write_corrs('{}{}_mm.csv'.format(i, j), corr)

        # do shear shear correlations
        # cross correlations are symmetric
        for (keyi, groupi), i in zip(self.source_groups, range(0, 3)):
            for (keyj, groupj), j in zip(self.source_groups, range(0, 3)):
                if (j > i):
                    continue
                # catch the auto correlation
                elif i == j:
                    corr = self.shearshear(groupi)
                    self.io.write_corrs('{}{}_gg.csv'.format(i, j), corr)
                # cross correlations
                else:
                    corr = self.shearshear(groupi, cat2=groupj)
                    self.io.write_corrs('{}{}_gg.csv'.format(i, j), corr)
    
        # tangential shear correlations
        # not symmetric
        for (keyl, groupl), i in zip(self.lens_groups, range(0, 3)):
            for (keys, groups), j in zip(self.source_groups, range(0, 3)):
                corr = self.gammat(groupl, groups)
                self.io.write_corrs('{}{}_gm.csv'.format(i, j), corr)
        return

    def wtheta(self, table, table2=None):
        '''calculate position position correlation'''
        #setup correlation objects
        dd = treecorr.NNCorrelation(min_sep=1.0, max_sep=80, nbins=10, sep_units='arcmin', bin_slop=.01)
        rand = treecorr.Catalog(ra=self.randoms['ra'].values, dec=self.randoms['dec'].values,
                                ra_units='radians', dec_units='radians')
        rr = treecorr.NNCorrelation(min_sep=1.0, max_sep=80, nbins=10, sep_units='arcmin', bin_slop=.01)
        dr = treecorr.NNCorrelation(min_sep=1.0, max_sep=80, nbins=10, sep_units='arcmin', bin_slop=.01)

        #deal with second catalog if need be
        if table2 is not None:
            cat = self.io.df_to_corr(table)
            cat2 = self.io.df_to_corr(table2)

            rd = treecorr.NNCorrelation(min_sep=1.0, max_sep=80, nbins=10, sep_units='arcmin', bin_slop=.01)

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
            
    def pos_pos_corr(self, pos1,pos2,posr1,posr2,w1=None,w2=None,same_zshell=False,same_cell=False,unique_encounter=False,num_threads=0):
        nbins    = 10
        min_sep  = 0.025 # 3 arcmin
        max_sep  = 1.5 # 90 arcmin 
        bin_size = (max_sep-min_sep)/nbins # roughly
        bin_slop = 0.01 # 0.1 -> 0.05 # 2pt_pipeline for des used bin_slop: 0.01 here: https://github.com/des-science/2pt_pipeline/blob/master/pipeline/twopt_pipeline.yaml
        # num_threads = 5 #None #0 # should query the number of cpus your computer has
        logger = None

        if same_zshell and same_cell: # auto
            ra, dec = pos1 # either 1 or 2 works, they're the same
            ra_R, dec_R = posr1
            w = w1
            cat   = treecorr.Catalog(ra=ra, dec=dec, ra_units='degrees', dec_units='degrees', w=w) 
            cat_R = treecorr.Catalog(ra=ra_R, dec=dec_R, ra_units='radians', dec_units='radians')

        else: #cross
            ra1, dec1 = pos1
            ra2, dec2 = pos2
            ra_R1, dec_R1 = posr1
            ra_R2, dec_R2 = posr2
            cat1   = treecorr.Catalog(ra=ra1, dec=dec1, ra_units='degrees', dec_units='degrees', w=w1) 
            cat2   = treecorr.Catalog(ra=ra2, dec=dec2, ra_units='degrees', dec_units='degrees', w=w2)
            cat1_R = treecorr.Catalog(ra=ra_R1, dec=dec_R1, ra_units='radians', dec_units='radians')
            cat2_R = treecorr.Catalog(ra=ra_R2, dec=dec_R2, ra_units='radians', dec_units='radians')

        DD = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, bin_slop=bin_slop,
                                    sep_units='degrees', logger=logger)
        RR = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, bin_slop=bin_slop,
                                    sep_units='degrees', logger=logger)
        DR = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, bin_slop=bin_slop, 
                                    sep_units='degrees', logger=logger)
        RD = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins, bin_slop=bin_slop, 
                                    sep_units='degrees', logger=logger)

        if same_zshell and same_cell: # auto
            # same z, same pix
            DD.process_auto(cat,num_threads=num_threads)
            RR.process_auto(cat_R,num_threads=num_threads)
            DR.process_cross(cat, cat_R,num_threads=num_threads)
            RD = DR.copy() 
        elif same_zshell:
            if unique_encounter: # the following two counts shouldn't be doubled up cuz they're the same in both directions
                DD.process_cross(cat1, cat2,num_threads=num_threads)
                RR.process_cross(cat1_R, cat2_R,num_threads=num_threads)
            else:
                DR.process_cross(cat1, cat2_R,num_threads=num_threads)
                DR.process_cross(cat2, cat1_R,num_threads=num_threads)
                RD = DR.copy()
        else: # different z  (can have different/same pix) when 2 cats have diff zshells it is enough to make them different even within the same hp pix
            DD.process_cross(cat1, cat2,num_threads=num_threads) # metric='Rperp')
            RR.process_cross(cat1_R, cat2_R,num_threads=num_threads)
            DR.process_cross(cat1, cat2_R,num_threads=num_threads)
            RD.process_cross(cat1_R, cat2,num_threads=num_threads) 

        DD.finalize()
        RR.finalize()
        DR.finalize()
        RD.finalize()

        theta = np.exp(DD.meanlogr)
        wtheta, var = DD.calculateXi(RR,DR,RD) # var gives the Poisson error
        err_poisson = np.sqrt(var)

        return {"xi":wtheta, "sig":err_poisson, "r":theta}

    def gammat(self, lens, sources):
        '''calculate tangential shear correlation'''
        lens_corr = self.io.df_to_corr(lens, shears=True)
        source_corr = self.io.df_to_corr(sources, shears=True)
        rand = treecorr.Catalog(ra=self.randoms['ra'].values, dec=self.randoms['dec'].values,
                                ra_units='radians', dec_units='radians')
        #now make correlation functions
        GGL = treecorr.NGCorrelation(min_sep=0.1, max_sep=90, nbins=10, sep_units='arcmin')
        GGL.process(lens_corr, source_corr)

        # calculate random signal
        GGL_rand = treecorr.NGCorrelation(min_sep=0.1, max_sep=90, nbins=10, sep_units='arcmin')
        GGL_rand.process(rand, source_corr)
        return {'xi+': GGL.xi - GGL_rand.xi, 'xi-' : GGL.xi_im, 'r':np.exp(GGL.meanlogr), 'sig':np.sqrt(GGL.varxi)}

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
        sig = np.sqrt(gg.varxi)

        return {'xip': xip, 'xim': xim, 'r':r, 'sig':sig}
