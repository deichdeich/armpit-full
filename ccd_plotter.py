# This is for plotting the ccd that shows the dereddening method
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

isochrones = [isochrones[np.where(isochrones["M_ini"]<10)]][0]
iso_colors = np.rec.array([(isochrones["W336MAG"]-isochrones["F475MAG"]),
			   (isochrones["F475MAG"]-isochrones["F814MAG"])],
			   names=("Iso(336-475)","Iso(475-814)"))
iso_color_x = iso_colors["Iso(336-475)"]
iso_color_y = iso_colors["Iso(475-814)"]

#data
data = np.genfromtxt(filename, names=True, delimiter=",")

#z_data = (data[~np.isnan(data['Z'])])
z_data = data
z_data = z_data[z_data['F336WF475W_RAW']<-0.5]
z_data = z_data[z_data['AF475W']>-0.9]
print 'doing ccd...'

plt.figure()
plt.scatter(data['F336WF475W_RAW'],data['F475WF814W_RAW'],edgecolors='none', color='gray')
plt.scatter(z_data['F336WF475W_RAW'],z_data['F475WF814W_RAW'],edgecolors='none', c=z_data['AF475W'])
cb = plt.colorbar()


cb.set_clim([-0.5,3.5])
cb.set_label('A(F475W)')



plt.scatter(iso_color_x,iso_color_y,color='black')

plt.xlim(-2,0)
plt.ylim(-1,2)


plt.xlabel('F336W-F475W')
plt.ylabel('F475W-F814W')

plt.show()

#plt.savefig('ad_phat/plots/ccd_{}_{}.png'.format(month,day))

print 'ccd complete.'
