# This is for plotting the cmd that shows the metal fitting method
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import time


day = time.strftime("%d")
month = time.strftime("%m")
filename = 'ad_phat/metal_fit/full_data_metals_{}_{}.csv'.format(month,day)

## isochrone stuff (open, cut, arrange)
isochrones = fits.open("isochrones/new_padova_4Myr_finez.fits")
isochrones = isochrones[1].data

isochrones = [isochrones[np.where(np.logical_and(isochrones["M_ini"]<10,isochrones["M_ini"]>2))]][0]
iso_colors = np.rec.array([(isochrones["W336MAG"]-isochrones["F475MAG"]),
			   (isochrones["F475MAG"])],
			   names=("Iso(336-475)","Iso(475)"))
iso_color_x = iso_colors["Iso(336-475)"]
iso_mag_y = iso_colors["Iso(475)"]

#data
data = np.genfromtxt(filename, names=True, delimiter=",")
data = data[data['F336WF475W_RAW']<-0.5]

z_data = (data[~np.isnan(data['Z'])])

print 'doing cmd...'

plt.figure()
plt.scatter(data['F336WF475W_RAW'],data['F475W_RAW'],edgecolors='none', color='gray')
plt.scatter(z_data['F336WF475W_RAW'],z_data['F475W_RAW'],edgecolors='none', c=np.log10(z_data['Z']/0.019))
cb = plt.colorbar()
cb.set_clim(-1,3.5)
plt.scatter(iso_color_x,iso_color_y,s=1,color='red')

plt.xlim(-2,0)
plt.ylim(-6,0)

cb.set_label('log(Z/0.019)')
plt.xlabel('F336W-F475W')
plt.ylabel('F475W')

ax = plt.gca()
ax.invert_yaxis()

plt.savefig('ad_phat/plots/cmd_{}_{}.png'.format(month,day))
print 'cmd complete.'
