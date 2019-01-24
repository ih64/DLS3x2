from os.path import join
import numpy as np
import pandas as pd
import treecorr

class io():
    def __init__(self):
        self.catalog_path = '/global/homes/i/ihasan/shear_gp/data_access/catalogs/'
        self.out_path = '/global/homes/i/ihasan/shear_gp/data_access/3x2Pipe/buzzard_corrs/'

    def setup_lens(self, table):
        '''helper to read in lens input catalog
        '''
        mask = table['R'] < 21
        mask &= table['R'] > 18
        return table[mask]

    def setup_source(self, table):
        '''helper to read in source input catalog
        '''
        # magnitude cuts
        mask = table['R'] < 24.5
        mask &= table['R'] > 21
        # shape cuts
        mask &= table['de'] < .3
        mask &= table['b'] > .4
#        mask &= table['status'] == 1
        table = table[mask].copy()

        # shear calibration
#        m_gamma = 6*np.power(10.0 ,-4) * np.power(table['R'] - 20, 3.26) + 1.04
#        table.loc[:,'e1'] *= m_gamma
#        table.loc[:,'e2'] *= m_gamma
        return table

    def df_to_corr(self, table, shears=False):
        """
        turn an pd dataframe into a treecorr catalog
        anticipating certain format form astropy catalog
        shape catalogs should be weigheted by shape noise
        position catalogs should not
        """
        if shears:
            cat = treecorr.Catalog(ra=table['alpha'].values, dec=table['delta'].values,
                                   ra_units='deg', dec_units='deg', g1=table['e1'], g2=table['e2'],
                                   w=1/(.25**2 + table['de']**2))

        else:
            cat = treecorr.Catalog(ra=table['alpha'].values, dec=table['delta'].values,
                                   ra_units='deg', dec_units='deg')

        return cat


    def read_randoms(self, catalog_name):
        '''
        quick helper function to return the catalog of randoms
        '''
        ext = catalog_name.split('.')[-1]
        if ext == 'hdf':
            rand_cat = pd.read_hdf(join(self.catalog_path, catalog_name), '/data')
            return rand_cat
        elif ext == 'csv':
            rand_cat = pd.read_csv(join(self.catalog_path, catalog_name))
            return rand_cat
        else:
            print('unrecognized file extension, exiting')
            return

    def write_corrs(self, out_name, corr):
        '''
        quick function to write outputs of correlation functions
        '''
        df = pd.DataFrame(data=corr)
        output_path = join(self.out_path, out_name)
        df.to_csv(output_path)
        return
