from os.path import join
import pandas as pd
import treecorr

class io():
    def __init__(self):
        self.catalog_path = '/global/cscratch1/sd/ihasan/'
        self.randoms_path = '/global/homes/i/ihasan/shear_gp/data_access/catalogs'
        self.out_path = '/global/homes/i/ihasan/shear_gp/data_access/3x2Pipe/corrs'

    def setup_lens(self):
        '''helper to read in lens input catalog
        '''
        lens_input = pd.read_hdf(join(self.catalog_path, 'lens_input.hdf'), '/data')
        return lens_input

    def setup_source(self):
        '''helper to read in source input catalog
        '''
        source_input = pd.read_hdf(join(self.catalog_path, 'source_input.hdf'), '/data')
        return source_input

    def df_to_corr(self, table, weight_key, shears=False):
        """
        turn an pd dataframe into a treecorr catalog
        anticipating certain format form astropy catalog
        shape catalogs should be weigheted by shape noise
        position catalogs should not
        """
        if shears:
            cat = treecorr.Catalog(ra=table['alpha'].values, dec=table['delta'].values,
                                   ra_units='deg', dec_units='deg', g1=table['e1'], g2=table['e2'],
                                   w=table[weight_key]/(.25**2 + table['de']**2))

        else:
            cat = treecorr.Catalog(ra=table['alpha'].values, dec=table['delta'].values,
                                   ra_units='deg', dec_units='deg', w=table[weight_key])

        return cat


    def read_randoms(self, catalog_name):
        '''
        quick helper function to return the catalog of randoms
        '''
        rand_cat = pd.read_csv(join(self.randoms_path, catalog_name))
        return rand_cat

    def write_corrs(self, out_name, corr):
        '''
        quick function to write outputs of correlation functions
        '''
        df = pd.DataFrame(data=corr)
        output_path = join(self.out_path, out_name)
        df.to_csv(output_path)
        return
