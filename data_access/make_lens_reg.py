from astropy import units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord

prefix = ['F1','F2','F3','F4','F5']

for f in prefix:
    t = ascii.read('catalogs/'+f+'_all.csv')
    lens_mask = t['r'] < 22
    lens_mask &= t['z_b'] > .25
    lens_mask &= t['z_b'] < .8

    t = t[lens_mask]

    sc = SkyCoord(ra=t['alpha']*u.degree, dec=t['delta']*u.degree)

    with open('regfiles/'+f+'_lens.reg','w') as f:
        header = '# Region file format: DS9 version 4.1\n'
        header += 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
        header += 'fk5\n'
        f.write(header)
        for c in sc:
            cord_str = c.to_string('hmsdms', sep=':')
            ra, dec = cord_str.split()
            ds9_str = 'circle({},{},12.0")\n'.format(ra, dec)
            f.write(ds9_str)

