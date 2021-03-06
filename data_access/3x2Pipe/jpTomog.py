from os.path import exists
import numpy as np
import pandas as pd
import yaml
import pdb

class Tomography():
    def __init__(self, filename='tomog.yaml'):
        ''' do something '''
        # read in the file paths and such
        # TODO
        # right now you have to make sure the lens and source cuts
        # are the same in the sql query to make the lens/source_pz file
        # and in source/lens cuts on the phot_shape file. should add checks
        # to ensure that is the case
        with open(filename, 'r') as f:
            config = yaml.load(f)
            self.lens_pz_file = config['lens_pz_file']
            self.source_pz_file = config['source_pz_file']
            self.phot_shape_file = config['phot_shape_file']
            self.pz_sel_thresh = config['pz_sel_thresh']
            self.lens_dict = config['lens']
            self.source_dict = config['source']

        file_names = [self.lens_pz_file, self.source_pz_file, self.phot_shape_file]
        for file in file_names:
            assert exists(file), 'cannot find {} read in from config'.format(file)
            
    def run(self):
        ''' create lens and source catalogs you can use in 3x2 pt '''
        self.setup_catalog(self.phot_shape_file, self.lens_pz_file, 'lens')
        self.setup_catalog(self.phot_shape_file, self.source_pz_file, 'source')
        return
    
    def setup_catalog(self, phot_file, pz_file, cuts):
        ''' make selection cuts on phot_shape file and join to the 
        lens_pz file. use weights in the lens_pz file to place galaxies in bins
        '''
        phot = pd.read_hdf(phot_file, '/data')
        pz = pd.read_hdf(pz_file, '/data')

        # toss out galaxies that have low significance of being in any bins
        mask = pz['w1'] > self.pz_sel_thresh
        mask |= pz['w2'] > self.pz_sel_thresh
        mask |= pz['w3'] > self.pz_sel_thresh
        pz = pz[mask].copy()
        
        # make the photometry cuts
        phot = self.phot_cuts(phot, cuts)
        
        # de-dup the photometry catalog
        phot = self.dedup(phot)
        pz = self.dedup(pz)
        
        # join datafames on the objid
        joined_df = phot.set_index('objid').join(pz.set_index('objid'),
                                                 how='inner',lsuffix='_l', rsuffix='_r')
        # assign bins to galaxies
        bin1 = joined_df['w1'] > self.pz_sel_thresh
        bin2 = joined_df['w2'] > self.pz_sel_thresh
        bin3 = joined_df['w3'] > self.pz_sel_thresh

        joined_df.loc[bin1, 'z_bin'] = 0
        joined_df.loc[bin2, 'z_bin'] = 1
        joined_df.loc[bin3, 'z_bin'] = 2

        # toss out any galaxies that got assigned to no bins
        bin_mask = joined_df['z_bin'] == -1
        joined_df = joined_df[~bin_mask].copy()
        # save the catalog, trimming columns
        if cuts == 'lens':
            joined_df.loc[:, ['subfield_l','Bdered', 'Vdered', 'Rdered', 'zdered',
                              'alpha', 'delta', 'z_b', 'z_bin', 'w1', 'w2', 'w3']].to_hdf(self.lens_dict['lens_out'], '/data')
        elif cuts == 'source':
            joined_df.loc[:, ['subfield_l', 'Bdered', 'Vdered', 'Rdered', 'zdered',
                           'alpha', 'delta', 'z_b', 'z_bin', 'w1', 'w2', 'w3','e1', 'e2', 'de']].to_hdf(self.source_dict['source_out'], '/data')
        return

    def phot_cuts(self, table, cuts):
        ''' take in the photometry_source catalog and do either lens or source cuts'''
        assert cuts in ['lens', 'source'], \
            'cuts needs to be either "lens" or "source", you typed {}'.format(cuts)

        if cuts == 'lens':
            mask = table['Rdered'] <= self.lens_dict['mag_high']
            mask &= table['Rdered'] > self.lens_dict['mag_low']
            table = table.loc[mask].copy()

        elif cuts == 'source':
            # magnitude cuts
            mask = table['Rdered'] < self.source_dict['mag_high']
            mask &= table['Rdered'] > self.source_dict['mag_low']
            if self.source_dict['calibrate']:
                # shape cuts
                mask &= table['de'] < self.source_dict['de']
                mask &= table['b'] > self.source_dict['b']
                mask &= table['status'] == 1
                table = table.loc[mask].copy()

                # shear calibration
                m_gamma = 6*np.power(10.0, -4) * np.power(table['Rdered'] - 20, 3.26) + 1.04
                table.loc[:, 'e1'] *= m_gamma
                table.loc[:, 'e2'] *= m_gamma
            else:
                table = table[mask].copy()
            
        # prepare place holder for bins
        table['z_bin'] = -1*np.ones(len(table), np.int8)
        return table

    def dedup(self, table):
        '''remove duplicate rows in the photometry catalog'''
        unique_objids, uniq_idx, counts = np.unique(table['objid'], return_index=True,
                                                    return_counts=True)

        unique_table = table.iloc[uniq_idx].copy()
        return unique_table

if __name__ == '__main__':
    tom = Tomography()
    tom.run()
