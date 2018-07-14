from astropy.io import ascii
from astropy.table import vstack
import pickle
from probes import getGGL 

def doall():
    F1_all = ascii.read('F1_all.csv')
    F2_all = ascii.read('F2_all.csv')
    F3_all = ascii.read('F3_all.csv')
    F4_all = ascii.read('F4_all.csv')
    F5_all = ascii.read('F5_all.csv')
    F_all = vstack([F1_all, F2_all, F3_all, F4_all, F5_all])

    L1_mask = ((F_all['z.z_b'] > .25) & (F_all['z.z_b'] < .45) \
               & (F_all['p.r'] < 22))
#    L1point5_mask = ((F_all['z.z_b'] > .25) & (F_all['z.z_b'] < .45) & (F_all['p.r'] > 22))
    L2_mask = ((F_all['z.z_b'] > .45) & (F_all['z.z_b'] < .6) \
               & (F_all['p.r'] < 22))
    L3_mask = ((F_all['z.z_b'] > .6) & (F_all['z.z_b'] < .8) \
               & (F_all['p.r'] < 22) & (F_all['p.r'] > 19.5))

    S1_mask = ((F_all['z.z_b'] > .5) & (F_all['z.z_b'] < .7) \
               & (F_all['p.r'] > 21.5) & (F_all['s.de'] < .25) & (F_all['s.b'] > .3))
    S2_mask = ((F_all['z.z_b'] > .7) & (F_all['z.z_b'] < .85) \
               & (F_all['p.r'] > 21.5) & (F_all['s.de'] < .25) & (F_all['s.b'] > .3))
    S3_mask = ((F_all['z.z_b'] > .85) & (F_all['p.r'] > 21.5) \
               & (F_all['s.de'] < .25) & (F_all['s.b'] > .3))


    L1S1GGL = getGGL(F_all[L1_mask], F_all[S1_mask])
    pickle.dump(L1S1GGL, open('debugging/L1-S1_gammat.pkl', 'wb'))
    L2S2GGL = getGGL(F_all[L2_mask], F_all[S2_mask])
    pickle.dump(L2S2GGL, open('debugging/L2-S2_gammat.pkl', 'wb'))
    L3S3GGL = getGGL(F_all[L3_mask], F_all[S3_mask])
    pickle.dump(L3S3GGL, open('debugging/L3-S3_gammat.pkl', 'wb'))

if __name__ == '__main__':
    doall()
