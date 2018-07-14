import treecorr
import numpy as np
from astropy.io import ascii
from astropy.table import vstack
import pickle


def doall():
    F1_all = ascii.read('F1_all.csv')
    F2_all = ascii.read('F2_all.csv')
    F3_all = ascii.read('F3_all.csv')
    F4_all = ascii.read('F4_all.csv')
    F5_all = ascii.read('F5_all.csv')
    F_all = vstack([F1_all, F2_all, F3_all, F4_all, F5_all])

    L1_mask = ((F_all['z.z_b'] > .25) & (F_all['z.z_b'] < .45) & (F_all['p.r'] < 22))
    L1point5_mask = ((F_all['z.z_b'] > .25) & (F_all['z.z_b'] < .45) & (F_all['p.r'] > 22))
    L2_mask = ((F_all['z.z_b'] > .45) & (F_all['z.z_b'] < .6) & (F_all['p.r'] < 22))
    L3_mask = ((F_all['z.z_b'] > .6) & (F_all['z.z_b'] < .8) & (F_all['p.r'] < 22) & (F_all['p.r'] > 19.5))

    S1_mask = ((F_all['z.z_b'] > .5) & (F_all['z.z_b'] < .7) & (F_all['p.r'] > 21.5) & (F_all['s.de'] < .25) & (F_all['s.b'] > .3))
    S2_mask = ((F_all['z.z_b'] > .7) & (F_all['z.z_b'] < .85) & (F_all['p.r'] > 21.5) & (F_all['s.de'] < .25) & (F_all['s.b'] > .3))
    S3_mask = ((F_all['z.z_b'] > .85) & (F_all['p.r'] > 21.5) & (F_all['s.de'] < .25) & (F_all['s.b'] > .3))

    L1 = F_all[L1_mask]
    L1point5 = F_all[L1point5_mask]
    L2 = F_all[L2_mask]
    L3 = F_all[L3_mask]
    S1 = F_all[S1_mask]
    S2 = F_all[S2_mask]
    S3 = F_all[S3_mask]

    S1_mbias = 6*(10**-4) * np.sign(S1['p.r'] - 20) * np.abs(S1['p.r'] - 20)**3.26 + 1.036
    S2_mbias = 6*(10**-4) * np.sign(S2['p.r'] - 20) * np.abs(S2['p.r'] - 20)**3.26 + 1.036
    S3_mbias = 6*(10**-4) * np.sign(S3['p.r'] - 20) * np.abs(S3['p.r'] - 20)**3.26 + 1.036

    S1_tc = treecorr.Catalog(ra=S1['p.alpha'].data, dec=S1['p.delta'].data, w = S1['s.de'],
                         ra_units='deg', dec_units='deg', g1=S1['s.e1']*S1_mbias, g2=S1['s.e2']*S1_mbias)
    S2_tc = treecorr.Catalog(ra=S2['p.alpha'].data, dec=S2['p.delta'].data, w = S2['s.de'],
                         ra_units='deg', dec_units='deg', g1=S2['s.e1']*S2_mbias, g2=S2['s.e2']*S2_mbias)
    S3_tc = treecorr.Catalog(ra=S3['p.alpha'].data, dec=S3['p.delta'].data, w = S3['s.de'],
                         ra_units='deg', dec_units='deg', g1=S3['s.e1']*S3_mbias, g2=S3['s.e2']*S3_mbias)


 
    NGobj_kwargs = {"min_sep":0.1, "max_sep":90, 'nbins':10, 'sep_units':'arcmin'}
    

    L1S1r = []
    L1S1gammat = []
    L1S1gammax = []
    L1size = len(L1)
    for i in range(0,2000):
        resampled_idx = np.random.randint(0,L1size,L1size)
        L1_tc = treecorr.Catalog(ra=L1['p.alpha'].data[resampled_idx], dec=L1['p.delta'].data[resampled_idx], 
                         ra_units='deg', dec_units='deg')
        NGobj = treecorr.NGCorrelation(**NGobj_kwargs)
        NGobj.process(L1_tc, S1_tc)
        L1S1r.append(np.exp(NGobj.meanlogr))
        L1S1gammat.append(NGobj.xi)
        L1S1gammax.append(NGobj.xi_im)
        
    pickle.dump({'meanr':np.array(L1S1r), 'gammat':np.array(L1S1gammat), 'gammax':np.array(L1S1gammax)},
            open('new/L1-S1_gammat.pkl','wb'))
    print('finished L1S1')


    L1S2r = []
    L1S2gammat = []
    L1S2gammax = []
    for i in range(0,2000):
        resampled_idx = np.random.randint(0,L1size,L1size)
        L1_tc = treecorr.Catalog(ra=L1['p.alpha'].data[resampled_idx], dec=L1['p.delta'].data[resampled_idx], 
                         ra_units='deg', dec_units='deg')
        NGobj = treecorr.NGCorrelation(**NGobj_kwargs)
        NGobj.process(L1_tc, S2_tc)
        L1S2r.append(np.exp(NGobj.meanlogr))
        L1S2gammat.append(NGobj.xi)
        L1S2gammax.append(NGobj.xi_im)

    pickle.dump({'meanr':np.array(L1S2r), 'gammat':np.array(L1S2gammat), 'gammax':np.array(L1S2gammax)},
            open('new/L1-S2_gammat.pkl','wb'))

    print('finished L1 S2')

    L1S3r = []
    L1S3gammat = []
    L1S3gammax = []
    for i in range(0,2000):
        resampled_idx = np.random.randint(0,L1size,L1size)
        L1_tc = treecorr.Catalog(ra=L1['p.alpha'].data[resampled_idx], dec=L1['p.delta'].data[resampled_idx], 
                            ra_units='deg', dec_units='deg')
        NGobj = treecorr.NGCorrelation(**NGobj_kwargs)
        NGobj.process(L1_tc, S3_tc)
        L1S3r.append(np.exp(NGobj.meanlogr))
        L1S3gammat.append(NGobj.xi)
        L1S3gammax.append(NGobj.xi_im)

    pickle.dump({'meanr':np.array(L1S3r), 'gammat':np.array(L1S3gammat), 'gammax':np.array(L1S3gammax)},
            open('new/L1-S3_gammat.pkl','wb'))

    print('finished L1 S3')


    L1point5S1r = []
    L1point5S1gammat = []
    L1point5S1gammax = []
    L1point5size = len(L1point5)
    for i in range(0,2000):
        resampled_idx = np.random.randint(0,L1point5size,L1point5size)
        L1point5_tc = treecorr.Catalog(ra=L1point5['p.alpha'].data[resampled_idx], dec=L1point5['p.delta'].data[resampled_idx], 
                         ra_units='deg', dec_units='deg')
        NGobj = treecorr.NGCorrelation(**NGobj_kwargs)
        NGobj.process(L1point5_tc, S1_tc)
        L1point5S1r.append(np.exp(NGobj.meanlogr))
        L1point5S1gammat.append(NGobj.xi)
        L1point5S1gammax.append(NGobj.xi_im)
        
    pickle.dump({'meanr':np.array(L1point5S1r), 'gammat':np.array(L1point5S1gammat), 'gammax':np.array(L1point5S1gammax)},
            open('new/L15-S1_gammat.pkl','wb'))
    print('finished L15S1')


    L1point5S2r = []
    L1point5S2gammat = []
    L1point5S2gammax = []
    for i in range(0,2000):
        resampled_idx = np.random.randint(0,L1point5size,L1point5size)
        L1point5_tc = treecorr.Catalog(ra=L1point5['p.alpha'].data[resampled_idx], dec=L1point5['p.delta'].data[resampled_idx], 
                         ra_units='deg', dec_units='deg')
        NGobj = treecorr.NGCorrelation(**NGobj_kwargs)
        NGobj.process(L1point5_tc, S2_tc)
        L1point5S2r.append(np.exp(NGobj.meanlogr))
        L1point5S2gammat.append(NGobj.xi)
        L1point5S2gammax.append(NGobj.xi_im)

    pickle.dump({'meanr':np.array(L1point5S2r), 'gammat':np.array(L1point5S2gammat), 'gammax':np.array(L1point5S2gammax)},
            open('new/L15-S2_gammat.pkl','wb'))

    print('finished L1 S2')

    L1point5S3r = []
    L1point5S3gammat = []
    L1point5S3gammax = []
    for i in range(0,2000):
        resampled_idx = np.random.randint(0,L1point5size,L1point5size)
        L1point5_tc = treecorr.Catalog(ra=L1point5['p.alpha'].data[resampled_idx], dec=L1point5['p.delta'].data[resampled_idx], 
                            ra_units='deg', dec_units='deg')
        NGobj = treecorr.NGCorrelation(**NGobj_kwargs)
        NGobj.process(L1point5_tc, S3_tc)
        L1point5S3r.append(np.exp(NGobj.meanlogr))
        L1point5S3gammat.append(NGobj.xi)
        L1point5S3gammax.append(NGobj.xi_im)

    pickle.dump({'meanr':np.array(L1point5S3r), 'gammat':np.array(L1point5S3gammat), 'gammax':np.array(L1point5S3gammax)},
            open('new/L15-S3_gammat.pkl','wb'))

    print('finished L15 S3')

    L2S2r = []
    L2S2gammat = []
    L2S2gammax = []
    L2size = len(L2)
    for i in range(0,2000):
        resampled_idx = np.random.randint(0,L2size,L2size)
        L2_tc = treecorr.Catalog(ra=L2['p.alpha'].data[resampled_idx], dec=L2['p.delta'].data[resampled_idx], 
                            ra_units='deg', dec_units='deg')
        NGobj = treecorr.NGCorrelation(**NGobj_kwargs)
        NGobj.process(L2_tc, S2_tc)
        L2S2r.append(np.exp(NGobj.meanlogr))
        L2S2gammat.append(NGobj.xi)
        L2S2gammax.append(NGobj.xi_im)

    pickle.dump({'meanr':np.array(L2S2r), 'gammat':np.array(L2S2gammat), 'gammax':np.array(L2S2gammax)},
            open('new/L2-S2_gammat.pkl','wb'))

    print('finished L2 S2 ')

    L2size = len(L2)
    L2S3r = []
    L2S3gammat = []
    L2S3gammax = []
    for i in range(0,2000):
        resampled_idx = np.random.randint(0,L2size,L2size)
        L2_tc = treecorr.Catalog(ra=L2['p.alpha'].data[resampled_idx], dec=L2['p.delta'].data[resampled_idx], 
                                 ra_units='deg', dec_units='deg')
        NGobj = treecorr.NGCorrelation(**NGobj_kwargs)
        NGobj.process(L2_tc, S3_tc)
        L2S3r.append(np.exp(NGobj.meanlogr))
        L2S3gammat.append(NGobj.xi)
        L2S3gammax.append(NGobj.xi_im)
        
    print('finished L2 S3')
    pickle.dump({'meanr':np.array(L2S3r), 'gammat':np.array(L2S3gammat), 'gammax':np.array(L2S3gammax)},
            open('new/L2-S3_gammat.pkl','wb'))



    NGobj_kwargs = {"min_sep":0.1, "max_sep":90, 'nbins':8, 'sep_units':'arcmin'}

    L3S3r = []
    L3S3gammat = []
    L3S3gammax = []
    L3size = len(L3)
    for i in range(0,3000):
        resampled_idx = np.random.randint(0,L3size,L3size)
        L3_tc = treecorr.Catalog(ra=L3['p.alpha'].data[resampled_idx], dec=L3['p.delta'].data[resampled_idx], 
                                 ra_units='deg', dec_units='deg')
        NGobj = treecorr.NGCorrelation(**NGobj_kwargs)
        NGobj.process(L3_tc, S3_tc)
        L3S3r.append(np.exp(NGobj.meanlogr))
        L3S3gammat.append(NGobj.xi)
        L3S3gammax.append(NGobj.xi_im)

    pickle.dump({'meanr':L3S3r, 'gammat':L3S3gammat, 'gammax':np.array(L3S3gammax)}, open('new/L3-S3_gammat.pkl','wb'))

    print('finished L3 S3')


if __name__ == '__main__':
    doall()
