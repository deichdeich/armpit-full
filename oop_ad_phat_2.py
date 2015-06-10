# This is me trying to make the analysis object-oriented.
# I don't really know what I'm doing.  Here's what I have
# in mind:  Have a "ad_phat" class with methods for reddening
# and metal fitting.  These methods would basically just
# take the code from reddening.py and metal_fitter.py.
# I think that having it object orientied will also make
# handling the simulation output easier, as I can (in theory)
# just have a simulation object which I can pass immediately
# to the ad_phat class.

from astropy.io import fits
import numpy as np
import scipy.interpolate as interpolate
from scipy.optimize import fsolve


class ad_phat(object):

	def __init__(self,data_path,iso_age):
		self.iso_age = iso_age
		self.data_path = data_path
		self.iso_path = 'isochrones/{}Myr_finez.fits'.format(self.iso_age)

		if self.data_path.endswith('.csv'):
			self.raw_data = np.genfromtxt(self.data_path, names=True, delimiter=',')
		elif self.data_path.endswith('.fits'):
			self.raw_data = fits.open(self.data_path)
			self.raw_data = self.raw_data[1].data
			self.raw_data = self.raw_data[np.where(self.raw_data['F814W_ERR']<0.25)]
			self.raw_data = self.raw_data[np.where(self.raw_data['F475W_ERR']<0.25)]
				
		self.isochrones = fits.open(self.iso_path)
		self.isochrones = self.isochrones[1].data
		self.isochrones = [self.isochrones[np.where(np.logical_and(self.isochrones["M_ini"]<10,self.isochrones["M_ini"]>2))]][0]
		self.iso_colors = np.rec.array([(self.isochrones["W336MAG"]-self.isochrones["F475MAG"]),
					   (self.isochrones["F475MAG"]-self.isochrones["F814MAG"])],
						names=("Iso(336-475)","Iso(475-814)"))

			
		self.data_colors = np.rec.array([(self.raw_data['F336W_VEGA']-self.raw_data['F475W_VEGA']),
				    (self.raw_data['F475W_VEGA']-self.raw_data['F814W_VEGA'])],
				    names = ("Data(336-475)","Data(475-814)"))

		self.iso_color_x = sorted(self.iso_colors["Iso(336-475)"])
		self.iso_color_y = sorted(self.iso_colors["Iso(475-814)"])

		self.raw_x = self.data_colors["Data(336-475)"]
		self.raw_y = self.data_colors["Data(475-814)"]

		self.fit_coeffs = np.polyfit(self.iso_color_x,self.iso_color_y,3)

	def redvec(Rv_value):
		if Rv_value == 3.1:
			return (0.446,0.610)
		elif Rv_value == 5:
		        return (0.208,0.467)
		else:
		        return ValueError("{} is not a valid Rv value".format(Rv_value))


	Rv = 3.1
	E336m475,E475m814 = redvec(Rv)


	def IsochroneFit(self,x):
		return (self.fit_coeffs[0]*x**3 + self.fit_coeffs[1]*x**2 + self.fit_coeffs[2]*x + self.fit_coeffs[3])
	

	def extinction_func(self,x,x1,x2):
		return ((self.E475m814/self.E336m475)*(x-x1)+x2)

	def findIntersection(self,isochrone_func,extinction_func,x1,x2):
		return (lambda x: isochrone_func(x)-extinction_func(x,x1,x2))

	def calc_extinction(self,starnum):
		mag475 = self.raw_data['F475W_VEGA'][starnum]
		r1 = self.raw_x[starnum]
		r2 = self.raw_y[starnum]

		intersect_function = self.findIntersection(self.IsochroneFit,self.extinction_func,r1,r2)
		intersectpoint = fsolve(intersect_function,-1.5)
		

		newx = intersectpoint
		newy = self.IsochroneFit(intersectpoint)
		dx = r1 - newx
		dy = r2 - newy

		a475 = dx/self.E336m475
		new475 = mag475 - 24.4 - a475
		return(newx,newy,new475,a475)
			
