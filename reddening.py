import photometry_errors as perr
from astropy.io import fits
import numpy as np
import scipy.interpolate as interpolate
from scipy.optimize import fsolve


## isochrone stuff (open, cut, arrange)
isochrones = fits.open("isochrones/new_padova_4Myr_finez.fits")
isochrones = isochrones[1].data

isochrones = [isochrones[np.where(np.logical_and(isochrones["M_ini"]<10,isochrones["M_ini"]>2))]][0]

iso_colors = np.rec.array([(isochrones["W336MAG"]-isochrones["F475MAG"]),
	(isochrones["F475MAG"]-isochrones["F814MAG"])],
	names=("Iso(336-475)","Iso(475-814)"))

iso_color_x = iso_colors["Iso(336-475)"]
iso_color_y = iso_colors["Iso(475-814)"]
	
## phat stuff(open, arrange)
'''
phatdata_full = fits.open("raw_data/brightms_6filt.fits")
phatdata_raw = phatdata_full[1].data
phatdata_814cut = phatdata_raw[np.where(phatdata_raw["F814W_ERR"]<0.25)]
phatdata = phatdata_814cut[np.where(phatdata_814cut["F475W_ERR"]<0.25)]
'''
phatdata = np.genfromtxt('sim_06_01.csv',names=True,delimiter=',')


data_colors = np.rec.array([(phatdata["F336W_VEGA"]-phatdata["F475W_VEGA"]),
    (phatdata["F475W_VEGA"]-phatdata["F814W_VEGA"])],
	names=("Data(336-475)","Data(475-814)"))

# raw data colors
raw_x = data_colors["Data(336-475)"]
raw_y = data_colors["Data(475-814)"]


## These are the TRUE, CORRECT dereddening vectors
def redvec(Rv_value):
    if Rv_value == 3.1:
    	return (0.446,0.610)
    elif Rv_value == 5:
    	return (0.208,0.467)
    else:
    	return ValueError("{} is not a valid Rv value".format(Rv_value))

Rv = 3.1
E336m475,E475m814 = redvec(Rv)

iso_color_x = sorted(iso_color_x)
iso_color_y = sorted(iso_color_y)
fit_coeffs = np.polyfit(iso_color_x,iso_color_y,3)

def IsochroneFit(x):
    return (fit_coeffs[0] * x ** 3 + fit_coeffs[1] * x ** 2 + fit_coeffs[2] * x + fit_coeffs[3])


def extinction_func(x,x1,x2):
    return ((E475m814/E336m475)*(x-x1)+x2)


def findIntersection(isochrone_func,extinction_func,x1,x2):
    return (lambda x: isochrone_func(x)-extinction_func(x,x1,x2))

def mag_error(starnum):
    # error on the magnitudes
    err336,err475,err814 = perr.get_in_out(phatdata[starnum])

    # error on the color
    err336475 = np.sqrt(err336**2+err475**2)
    err475814 = np.sqrt(err475**2+err814**2)

    # distribution of possible colors
    dx = np.random.normal(0,err336475,100).tolist()
    dy = np.random.normal(0,err475814,100).tolist()
    deltas = zip(dx,dy)
    return deltas


def color_error(starnum):
    
    mag475 = phatdata["F475W_VEGA"][starnum]
    
    #deltas = mag_error(starnum)
    
    r1 = raw_x[starnum]
    r2 = raw_y[starnum]
    intersect = findIntersection(IsochroneFit,extinction_func,r1,r2)
    intersectpoint = fsolve(intersect,-1.5)
    newx = (intersectpoint)
    newy = (IsochroneFit(intersectpoint))
    dx = (r1-intersectpoint)
    dy = (r2-IsochroneFit(intersectpoint))
        
    a475 = dx/E336m475
    new475 = mag475-24.4-a475
    
    return (phatdata["RA"][starnum],
    phatdata["DEC"][starnum],
    r1,
    r2,
    mag475,
    newx[0],
    newy[0],
    new475[0],
    a475[0])
    
import time
month = time.strftime("%m")
day = time.strftime("%d")
filename = "ad_phat/dereddened_populations/full_data_{}_{}.csv".format(month,day)
with open(filename,"wb") as outfile:
    outfile.write("RA,DEC,F336W-F475W_RAW,F475W-F814W_RAW,F475W_RAW,F336W-F475W,F475W-F814W,F475W,A(F475W)\n")


def write_to_file(idx):
    with open(filename,"ab") as outfile:
        data = color_error(idx)
        outfile.write("{}\n".format(','.join(map(str,(data)))))
        if idx%1000 == 0:
            print "{} written".format(idx)

def mp_write():
    import multiprocessing as mp
    starlist = xrange(len(phatdata))
    num_proc = 4
    pool = mp.Pool(num_proc)
    pool.map(write_to_file,starlist)
