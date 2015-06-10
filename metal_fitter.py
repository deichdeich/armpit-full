import numpy as np
from astropy.io import fits
import scipy.interpolate as intp
import time


day = time.strftime("%d")
month = time.strftime("%m")
filename = "ad_phat/dereddened_populations/full_data_{}_{}.csv".format(month,day)

data = np.genfromtxt(filename, names=True, delimiter=",")

isochrones = fits.open("isochrones/4Myr_finez.fits")
isochrones = isochrones[1].data
isochrones = isochrones[isochrones["M_ini"]<20]

isocolor = isochrones["W336MAG"] - isochrones["F475MAG"]
isomag = isochrones["F475MAG"]
points = zip(isocolor,isomag)
zvalues = isochrones["Z"]

metals = intp.LinearNDInterpolator(points,zvalues)


datacolor = data["F336WF475W"]
#datacolorsig = data["F336WF475W_sigma"]
datamag = data["F475W"]
#datamagsig = data["F475W_sigma"]


def get_dxdy(idx):
    dx = np.random.normal(0,datacolorsig[idx],100).tolist()
    dy = np.random.normal(0,datamagsig[idx],100).tolist()
    deltas = zip(dx,dy)
    return deltas
    
def get_z(idx):
    color = datacolor[idx]
    mag = datamag[idx]
    z_no_error = metals(color,mag)
    return z_no_error



def write_to_file():
    with open("ad_phat/metal_fit/full_data_metals_{}_{}.csv".format(month,day),"wb") as outfile:
        outfile.write("RA,DEC,F336W-F475W_RAW,F475W-F814W_RAW,F475W_RAW,F336W-F475W,F475W-F814W,F475W,A(F475W),Z\n")
    with open("ad_phat/metal_fit/full_data_metals_{}_{}.csv".format(month,day),"ab") as outfile:
        for idx in xrange(len(data)):
	    print idx
            ra = data["RA"][idx]
            dec = data["DEC"][idx]
	    r1 = data["F336WF475W_RAW"][idx]
	    r2 = data["F475WF814W_RAW"][idx]
	    f475 = data["F475W_RAW"][idx]
	    r1new = data["F336WF475W"][idx]
	    r2new = data["F475WF814W"][idx]
	    f475new = data["F475W"][idx]
	    a475 = data["AF475W"][idx]
            zinfo = get_z(idx)
            outfile.write("{},{},{},{},{},{},{},{},{},{}\n".format(ra,dec,r1,r2,f475,r1new,r2new,f475new,a475,zinfo))
            if idx % 1000 == 0:
                print "{} written".format(idx)    
