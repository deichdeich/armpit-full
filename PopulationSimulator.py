import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.interpolate as intp
import os

'''
random_inf is a  funny little function; probably maxes out my personal record
for ratio of time spent thinking to number of lines written.

Input: int, number of stars you want sampled from the IMF
Returns: array of mass of each star determined by a dN/dM = M^(-2.35)
salpeter 1955 IMF.

intconst comes from integrating (1-2.35 = -1.35)
normconst comes from normalizing for the desired range (Integral from 1 to 100
equals 1.)
'''
def random_imf(size):
    x = np.random.random(size)
    intconst = 1.35
    normconst = 1.3527
    masses = (intconst/normconst)*(x**(1/(-intconst)))
    masses = masses[masses < 100]
    return masses
    
def get_isochrone(age):
    isochrone = fits.open('isochrones/{}Myr_finez.fits'.format(age))
    isochrone = isochrone[1].data
    #isochrone = [isochrone[np.where(np.logical_and(isochrone["M_ini"]<20,isochrone["M_ini"]>0))]][0]
    return isochrone

def interp(field,age):
    isochrone = get_isochrone(age)
    intp_function = intp.interp1d(isochrone['M_INI'],
        isochrone['{}'.format(field)],
        bounds_error=False)
    return intp_function

def add_grid(how_old,data):
    '''
    In color-magnitude space, add_grid adds a range of av's to the F475W
    magnitude.  Then, using the coefficients in extinction_coeffs.dat, it figures out
    the color that corresponds to that new magnitude in both the W336-F475W and the
    F475W-F814W colors.  Then it gets the corresponding magnitudes in the W336 and
    F814W filters.
    
    It does this by iterating over the list of av.  For each iteration, it calculates
    the relevant magnitudes, compiles them all in an array, and then appends that
    array onto the data file.
    '''

    # the extinction vectors from extinction_coeffs.dat
    redvecs = {'w336f475':[(-1)*(1.225/0.446)], 'f475f814':[(-1)*(1.225/0.61)]}
    
    # I don't really need to declare these as new variables, but it helps me keep track.
    new_av = data['A(F475W)']
    new_f475 = data['F475MAG']
    
    # I do need to declare this, though, so that each iteration adds to the one before it,
    # and doesn't just overwrite
    avgrid = data
    
    
    # the range of Av's
    f475_av_range = np.arange(0, 3.5, 0.025)
    for av in f475_av_range:
        new_f475 = data['F475MAG'] + av
        new_av = data['A(F475W)'] + av
        
        # calculated with (y2-y1) = m*(x2-x1) => x2 = [(y2-y1) + m*x1]/m
        new_w336f475 = (((-new_f475+data['F475MAG']) + (redvecs['w336f475'] * data['W336MAG-F475MAG']))/redvecs['w336f475'])
        new_f475f814 = (((-new_f475+data['F475MAG']) + (redvecs['f475f814'] * data['F475MAG-F814MAG']))/redvecs['f475f814'])
        new_w336 = new_w336f475 - new_f475
        new_f814 = -1 * (new_f475f814-new_f475)
        
        null_column = np.zeros(len(data))
        
        new_data = np.array(zip(null_column,
				null_column,
				data['MASS'],
                                new_w336,
                                new_f475,
                                new_f814,
                                new_w336f475,
                                new_f475f814,
                                new_av), dtype=[('RA',float),
						('DEC',float),
						('MASS',float),
                                                ('W336MAG',float),
                                                ('F475MAG',float),
                                                ('F814MAG',float),
                                                ('W336MAG-F475MAG',float),
                                                ('F475MAG-F814MAG',float),
                                                ('A(F475W)',float)])
        

        avgrid = np.append(avgrid,new_data)
    return avgrid



def add_real_distribution(how_old,data):
    '''
    Add Av from a distribution gotten from the original analysis.  This distribution
    has been selected from a green/blue region in the Dalcanton comparison map.
    '''
    redvecs = {'w336f475':[(-1)*(1.225/0.446)], 'f475f814':[(-1)*(1.225/0.61)]}
    
    av_dist_data = np.genfromtxt('good_av_distribution.csv',names=True,delimiter=',')
    av_dist_list = -av_dist_data['AF475W']
    
    sim_av = data['A(F475W)']
    sampled_av = np.random.choice(av_dist_list,len(sim_av))
    sim_av = sampled_av
    
    new_f475=data['F475MAG'] + sim_av
    
    # calculated with (y2-y1) = m*(x2-x1) => x2 = [(y2-y1) + m*x1]/m
    new_w336f475 = (((-new_f475+data['F475MAG']) + (redvecs['w336f475'] * data['W336MAG-F475MAG']))/redvecs['w336f475'])
    new_f475f814 = (((-new_f475+data['F475MAG']) + (redvecs['f475f814'] * data['F475MAG-F814MAG']))/redvecs['f475f814'])
    new_w336 = new_w336f475 - new_f475
    new_f814 = -1 * (new_f475f814-new_f475)
    null_column = np.zeros(len(sim_av))
    index_column = np.arange(len(sim_av))
    av_dist = np.array(zip(index_column,
                                index_column,
                                data['MASS'],
                                new_w336,
                                new_f475,
                                new_f814,
                                new_w336f475,
                                new_f475f814,
                                sim_av,
                                null_column,
                                null_column,
                                null_column,
                                null_column), dtype=[('RA',float),
                                                ('DEC',float),
                                                ('MASS',float),
                                                ('W336MAG',float),
                                                ('F475MAG',float),
                                                ('F814MAG',float),
                                                ('W336MAG-F475MAG',float),
                                                ('F475MAG-F814MAG',float),
                                                ('A(F475W)',float),
                                                ('F336WF475W_sigma',int),
                                                ('F475W_sigma',int),
                                                ('F475W_sigma2',int),
                                                ('F475W_sigma3',int)])
    return av_dist

def add_extinction(how_old,data, add_type='avgrid'):
    if add_type == 'avgrid':
        return add_grid(how_old,data)
    elif add_type == 'realav':
        return add_real_distribution(how_old,data)
    else:
        raise ValueError('{} not an option'.format(add_type))

def write_to_file(numstars,how_old,addtype='avgrid'):
    masslist = random_imf(numstars)
    ## the interpolation is meaningless for mass > 50
    masslist = masslist[masslist <= 50]
    w336_mag = interp('W336MAG', how_old)(masslist)
    f475_mag = interp('F475MAG', how_old)(masslist)
    f814_mag = interp('F814MAG', how_old)(masslist)
    w336f475 = w336_mag - f475_mag
    f475f814 = f475_mag - f814_mag
    
    # the original data points have no extinction
    av = np.zeros(numstars)
    no_av_data = np.array(zip(np.zeros(numstars),
    			np.zeros(numstars),
    			masslist,
                        w336_mag,
                        f475_mag,
                        f814_mag,
                        w336f475,
                        f475f814,
                        av),
                        dtype=[('RA',float),
			    ('DEC',float),
			    ('MASS',float),
                            ('W336MAG',float),
                            ('F475MAG',float),
                            ('F814MAG',float),
                            ('W336MAG-F475MAG',float),
                            ('F475MAG-F814MAG',float),
                            ('A(F475W)',float)])
    
    data_with_av = add_extinction(how_old,no_av_data,add_type=addtype)
    import time
    month = time.strftime('%m')
    day = time.strftime('%d')
    
        
    filename = 'simulation/{}myr/input_populations/{}_{}.csv'.format(how_old,month,day)
    np.savetxt(filename,
        data_with_av,
        fmt = '%.5f',
        delimiter = ',',
        comments = '',
        header = 'RA,DEC,MASS,F336W_VEGA,F475W_VEGA,F814W_VEGA,W336MAG-F475MAG,F475MAG-F814MAG,A(F475W)')
        
def simulate(nmbr,agespan=[4,30,100],addtype=['avgrid','realav']):
    if agespan != [4,30,100]:
        agespan = [agespan]
    if addtype != ['avgrid','realav']:
        addtype = [addtype]
    for age in agespan:
        for avtype in addtype:
            write_to_file(nmbr,age,addtype=avtype)
            print '{}Myr, {} finished'.format(age,avtype)
