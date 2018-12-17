import pickle
from astropy.io import ascii
from astropy.table import vstack
from probes import getGGL 
import numpy as np

def trim(f_all):
    f_group_by = f_all.group_by('subfield')
    f_groups = f_group_by.groups
    f_cleaned = []
    for sf in f_groups:
        ra_max = sf['alpha'].max()
        ra_min = sf['alpha'].min()
        delta_max = sf['delta'].max()
        delta_min = sf['delta'].min()

        ds = 5/60
        dra = ds/np.cos((np.pi/180.)*np.mean([delta_max, delta_min]))

        ra_max -= dra
        ra_min += dra
        delta_max -= ds
        delta_min += ds
    
        mask = sf['alpha'] < ra_max
        mask &= sf['alpha'] > ra_min
        mask &= sf['delta'] < delta_max
        mask &= sf['delta'] > delta_min
    
        f_cleaned.append(sf[mask])

    return vstack(f_cleaned)

def doall():
    F1_all = ascii.read('catalogs/F1_all.csv')
#    F1_all = trim(F1_all)

    F2_all = ascii.read('catalogs/F2_all.csv')
#    F2_all = trim(F2_all)

    F3_all = ascii.read('catalogs/F3_all.csv')
#    F3_all = trim(F3_all)

    F4_all = ascii.read('catalogs/F4_all.csv')
#    F4_all = trim(F4_all)

#    stars = ascii.read('catalogs/stars.csv')

    F5_all = ascii.read('catalogs/F5_all.csv')
#    F5_all = trim(F5_all)

    F_all = vstack([F1_all, F2_all, F3_all, F4_all, F5_all])

    L1_mask = ((F_all['z_b'] > .25) & (F_all['z_b'] < .45) \
               & (F_all['r'] < 21))
#    L1point5_mask = ((F_all['z.z_b'] > .25) & (F_all['z.z_b'] < .45) & (F_all['p.r'] > 22))
    L2_mask = ((F_all['z_b'] > .45) & (F_all['z_b'] < .6) \
               & (F_all['r'] < 21))
    L3_mask = ((F_all['z_b'] > .6) & (F_all['z_b'] < .8) \
               & (F_all['r'] < 21) & (F_all['r'] > 19.5))

    shape_mask = F_all['de'] <= .3
    shape_mask &= F_all['b'] > .6
#    shape_mask &= F_all['b'] < 3
    shape_mask &= np.sqrt(F_all['e2']**2 + F_all['e1']**2) < .6
#    shape_mask &= np.sqrt( 1 - np.square(F_all['b']) / np.square(F_all['a'])) < .9

    S1_mask = ((F_all['z_b'] > .5) & (F_all['z_b'] < .7) \
               & (F_all['r'] > 21) & (shape_mask))
    S2_mask = ((F_all['z_b'] > .7) & (F_all['z_b'] < .85) \
               & (F_all['r'] > 21) & (shape_mask))
    S3_mask = ((F_all['z_b'] > .85) & (F_all['r'] > 21) \
               & (shape_mask))

    kwargs = {'n_resample':1000, 'swap_test':False}
    L1S1GGL = getGGL(F_all[L1_mask], F_all[S1_mask], **kwargs)
    pickle.dump(L1S1GGL, open('draconian/L1-S1_gammat.pkl', 'wb'))
    L1S2GGL = getGGL(F_all[L1_mask], F_all[S2_mask], **kwargs)
    pickle.dump(L1S2GGL, open('draconian/L1-S2_gammat.pkl', 'wb'))
    L1S3GGL = getGGL(F_all[L1_mask], F_all[S3_mask], **kwargs)
    pickle.dump(L1S3GGL, open('draconian/L1-S3_gammat.pkl', 'wb'))

    L2S2GGL = getGGL(F_all[L2_mask], F_all[S2_mask], **kwargs)
    pickle.dump(L2S2GGL, open('draconian/L2-S2_gammat.pkl', 'wb'))
    L2S3GGL = getGGL(F_all[L2_mask], F_all[S3_mask], **kwargs)
    pickle.dump(L2S3GGL, open('draconian/L2-S3_gammat.pkl', 'wb'))

    L3S3GGL = getGGL(F_all[L3_mask], F_all[S3_mask], **kwargs)
    pickle.dump(L3S3GGL, open('draconian/L3-S3_gammat.pkl', 'wb'))


if __name__ == '__main__':
    doall()
