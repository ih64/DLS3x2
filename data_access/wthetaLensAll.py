from astropy.io import ascii
from probes import *
import pickle

def doall():
    L1 = {}
    L2 = {}
    L3 = {}

    for field in ['F1','F2','F3','F4','F5']:

        table = ascii.read(field+'_all.csv')

        L1_mask = ((table['z.z_b'] > .3) & (table['z.z_b'] < .45) & (table['p.r'] < 22))
        L2_mask = ((table['z.z_b'] > .45) & (table['z.z_b'] < .6) & (table['p.r'] < 22))
        L3_mask = ((table['z.z_b'] > .6) & (table['z.z_b'] < .8) & (table['p.r'] < 22) & (table['p.r'] > 19.5))
        [table.rename_column(c, c[2:]) for c in table.colnames]

        probe_list = calcProbes(table[L1_mask], field)
        L1[field] = probe_list
        
        probe_list = calcProbes(table[L2_mask], field)
        L2[field] = probe_list

        probe_list = calcProbes(table[L3_mask], field)
        L3[field] = probe_list

        print("finished field %s" % field)

    pickle.dump(L1,open('L1_wtheta.pkl','wb'))
    pickle.dump(L2,open('L2_wtheta.pkl','wb'))
    pickle.dump(L3,open('L3_wtheta.pkl','wb'))

    return

if __name__ == '__main__':
    doall()
