import numpy as np
import pylab as plt
from glob import glob
from astropy.io import fits

filepath = glob("../data/artificial_stars/SixFilt/*B09*F02.gst.*")[0]
astdata = fits.open(filepath)[1].data

#enable this to test the single-star reddening_with_errors.py
"""
phatdata = fits.open("brightms_6filt.fits")
phatdata = phatdata[1].data
"""

# xbins are global because you want to bin your sample data point as well
xbins = np.linspace(0,0.25,20)
dx = xbins[1]-xbins[0]
xcenters = xbins+dx

def binning(filt):
    if filt == "f336":
        m_in = astdata["F336W_IN"]
        m_out = astdata["F336W_VEGA"]
        err = astdata["F336W_ERR"]
    elif filt == "f475":
        m_in = astdata["F475W_IN"]
        m_out = astdata["F475W_VEGA"]
        err = astdata["F475W_ERR"]
    elif filt == "f814":
        m_in = astdata["F814W_IN"]
        m_out = astdata["F814W_VEGA"]
        err = astdata["F814W_ERR"]
    else:
        raise ValueError("{} is not an available filter".format(filt))
    
    m_out[m_out>90] = np.nan 
    err[err>90] = np.nan 

    # bias and compare
    bias = np.abs(m_in-m_out)
    comp = bias/err

    # use numpy
    idx = np.digitize(err,xbins)
    idx[idx==len(xbins)] = 0
    idx -= 1

    # group
    groups = {}
    for i,k in enumerate(idx):
        groups.setdefault(k,[]).append(i)
    
    return groups, bias, comp

"""
you're supposed to give it a DOLPHOT error of a data point and it will
give you the median of the corresponding bin of the in-out error.  You
can specify to give the mean instead of the median, and the comparison
value instead of the bias.
"""
def get_median(err, filt, func = np.median):
    index = np.digitize([err],xbins)
    index -= 1
    groups, bias, comp = binning(filt)
    group = groups[index[0]]
    in_out = func(bias[group])
    return in_out
    
    
def get_in_out(datapoint):
    
    err336 = datapoint["F336W_ERR"]
    err475 = datapoint["F475W_ERR"]
    err814 = datapoint["F814W_ERR"]
    
    # enable this to test the single-star reddening_with_errors.py
    """
    err336 = phatdata["F336W_ERR"][datapoint]
    err475 = phatdata["F475W_ERR"][datapoint]
    err814 = phatdata["F814W_ERR"][datapoint]
    """
    inout336 = get_median(err336,"f336")
    inout475 = get_median(err475,"f475")
    inout814 = get_median(err814,"f814")
    #print "DOLPHOT errors: {}(336), {}(475), {}(814)".format(err336, err475, err814)
    #print "In-Out errors: {}(336), {}(475), {}(814)".format(inout336, inout475, inout814)
    return inout336, inout475, inout814


    
