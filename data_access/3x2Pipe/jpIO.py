from os.path import join, exists
import pandas as pd
import yaml
import treecorr

class io():
    def __init__(self, filename='config.yaml'):
        # TODO make the config dict the attribute and use it later in the program
        # sanatize the yaml file
        config_keys = ['catalog_path', 'out_path', 'random_path', 'random_prefix']
        with open(filename,'r') as f:
            path_dict = yaml.load(f)
            for key in path_dict.keys():
                assert key in config_keys, 'key must be in {}'.format(config_keys)
                if 'path' in key:
                    assert exists(path_dict[key]), 'path {} does not exist'.format(path_dict[key])
                
        self.catalog_path = path_dict['catalog_path']
        self.out_path = path_dict['out_path']
        self.random_path = path_dict['random_path']
        self.random_prefix = path_dict['random_prefix']
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
            rand_cat = pd.read_hdf(join(self.random_path, catalog_name), '/data')
            rand_corr = treecorr.Catalog(ra=rand_cat['ra'].values, dec=rand_cat['dec'].values,
                                         ra_units='radians', dec_units='radians')
            return rand_corr
        
        elif ext == 'csv':
            rand_cat = pd.read_csv(join(self.random_path, catalog_name))
            rand_corr = treecorr.Catalog(ra=rand_cat['ra'].values, dec=rand_cat['dec'].values,
                                    ra_units='radians', dec_units='radians')
            return rand_corr
        
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
