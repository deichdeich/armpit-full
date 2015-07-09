import ad_phat_tools as apt
import numpy as np
from astropy.io import fits


artificial_star = apt.ArtificialInterpolatedMagnitudes()
masses = artificial_star.masslist(1500000)
print 'got masses for {} stars'.format(len(masses))

simage = 'big'

ra = np.zeros(len(masses))
dec = np.zeros(len(masses))
#age_arr = np.zeros(len(masses))+simage
age_arr = (100-4)*np.random.random(len(masses))+4
#z_arr = np.zeros(len(masses))+0.019
z_arr = (0.03-0.0001)*np.random.random(len(masses))+0.0001
f336w_arr = np.zeros(len(masses))
f475w_arr = np.zeros(len(masses))
f814w_arr = np.zeros(len(masses))

counter = 0
for starmass in masses:

    mags = artificial_star(age = age_arr[counter],mass = starmass,z_val = z_arr[counter],av = 0)
    f336w_arr[counter] = mags[0]
    f475w_arr[counter] = mags[1]
    f814w_arr[counter] = mags[2]
    ra[counter] = counter
    dec[counter] = counter
    print counter,' done'
    counter+=1
    
f336wf475w_arr = f336w_arr-f475w_arr
f475wf814w_arr = f475w_arr-f814w_arr

sim_data = fits.BinTableHDU.from_columns(
            [fits.Column(name='RA',format='E',array = ra),
             fits.Column(name='DEC',format = 'E',array = dec),
             fits.Column(name='MASS',format = 'E',array = masses),
             fits.Column(name='AGE',format = 'E',array = age_arr),
             fits.Column(name='Z',format = 'E',array = z_arr),
             fits.Column(name='F336W',format = 'E',array = f336w_arr),
             fits.Column(name='F475W',format = 'E',array = f475w_arr),
             fits.Column(name='F814W',format = 'E',array = f475w_arr),
             fits.Column(name='F336W-F475W',format = 'E',array = f336wf475w_arr),
             fits.Column(name='F475W-F814W',format = 'E',array = f475wf814w_arr)])

prihdr = fits.Header()
prihdr['IS_SIM'] = True
prihdr['SIM_AGE'] = simage
prihdu = fits.PrimaryHDU(header=prihdr)
hdulist = fits.HDUList([prihdu,sim_data])
hdulist.writeto('many_age_{}Myr.fits'.format(simage),clobber=True)
