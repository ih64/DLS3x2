import pickle
from astropy.io import ascii
from astropy.table import vstack
from probes import calcProbes

def doall():

    F1_all = ascii.read('catalogs/F1_all.csv')
    F2_all = ascii.read('catalogs/F2_all.csv')
    F3_all = ascii.read('catalogs/F3_all.csv')
    F4_all = ascii.read('catalogs/F4_all.csv')
    F5_all = ascii.read('catalogs/F5_all.csv')
    F_all = vstack([F1_all, F2_all, F3_all, F4_all, F5_all])
    
    L15_mask = ((F_all['z.z_b'] > .25) & (F_all['z.z_b'] < .45) & (F_all['p.r'] > 22)) 
#    L1_mask = ((F_all['z.z_b'] > .25) & (F_all['z.z_b'] < .45) & (F_all['p.r'] < 22))
#    L2_mask = ((F_all['z.z_b'] > .45) & (F_all['z.z_b'] < .6) & (F_all['p.r'] < 22))
#    L3_mask = ((F_all['z.z_b'] > .6) & (F_all['z.z_b'] < .8) & (F_all['p.r'] < 22) & \
#               (F_all['p.r'] > 19.5))


    #calculate wtheta
#    L1 = calcProbes(F_all[L1_mask], n_resample=1000)
#    pickle.dump(L1,open('debugging/L1_wtheta.pkl','wb'))
#    L2 = calcProbes(F_all[L2_mask], n_resample=1000)
#    pickle.dump(L2, open('debugging/L2_wtheta.pkl', 'wb'))
#    L3 = calcProbes(F_all[L3_mask], n_resample=1000)
#    pickle.dump(L3, open('debugging/L3_wtheta.pkl', 'wb'))
    L15 = calcProbes(F_all[L15_mask], n_resample=100)
    pickle.dump(L15, open('debugging/L1.5_wtheta.pkl','wb'))

    return

if __name__ == '__main__':
    doall()
