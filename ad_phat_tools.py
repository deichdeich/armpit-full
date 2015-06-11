from astropy.io import fits
import numpy as np
import scipy.interpolate as interpolate
from scipy.optimize import fsolve
import os.path
import time


class armpit(object):

    def __init__(self,data_path,iso_age):
        self.iso_age = iso_age
        self.iso_path = 'isochrones/{}Myr_finez.fits'.format(self.iso_age)
        self.isochrones = fits.open(self.iso_path)
        self.isochrones = self.isochrones[1].data
        self.isochrones = [self.isochrones[np.where(np.logical_and(self.isochrones["M_ini"]<10,
                                                                   self.isochrones["M_ini"]>2))]][0]
                                                                   
        self.iso_colors = np.rec.array([(self.isochrones["W336MAG"]-self.isochrones["F475MAG"]),
            (self.isochrones["F475MAG"]-self.isochrones["F814MAG"])],
            names=("Iso(336-475)","Iso(475-814)"))
                                        
        self.iso_color_x = sorted(self.iso_colors["Iso(336-475)"])
        self.iso_color_y = sorted(self.iso_colors["Iso(475-814)"])
        
        self.iso_mag = self.isochrones['F475MAG']
        self.zvalues = self.isochrones['Z']
        
        self.interp_points = zip(self.iso_color_x,self.iso_mag)
        self.metal_interp_function = interpolate.LinearNDInterpolator(self.interp_points,
                                                                      self.zvalues)
        
        
        self.data_path = data_path
        
        # get the path of the fits file    
        if type(self.data_path) is simulation:
            self.raw_file = self.data_path.filename
        elif self.data_path.endswith('.fits'):
            self.raw_file = self.data_path
        else:
            raise IOError('Data input must either be a fits file with containing the relevant columns or a simulation object.')
        
        self.raw_fits = fits.open(self.raw_file)
        
        self.data_header = self.raw_fits[0].header
        
        # check if it's a simulation output-- I *could* just do this when I do
        # if type(...) above, but I can conceive of having simulation output
        # that's not just an object.
        if 'IS_SIM' in self.data_header:
            self.is_sim = self.data_header['IS_SIM']
        else:
            self.is_sim = False
        
        # trim the data
        self.raw_data = self.raw_fits[1].data
        self.raw_data = self.raw_data[np.where(self.raw_data['F814W_ERR']<0.25)]
        self.raw_data = self.raw_data[np.where(self.raw_data['F475W_ERR']<0.25)]

            
        #self.data_colors = np.rec.array([(self.raw_data['F336W_VEGA']-self.raw_data['F475W_VEGA']),
                                        #(self.raw_data['F475W_VEGA']-self.raw_data['F814W_VEGA'])],
                                        #names = ("Data(336-475)","Data(475-814)"))

        

        #self.raw_x = self.data_colors["Data(336-475)"]
        #self.raw_y = self.data_colors["Data(475-814)"]

        self.fit_coeffs = np.polyfit(self.iso_color_x,self.iso_color_y,3)

        self.month = time.strftime('%m')
        self.day = time.strftime('%d')
        self.filename = 'armpit_out/armpit_out_{}_{}.csv'.format(self.month,self.day)
        
    '''
    Arguments: which columns you want, list
    Returns: data from inputted columns
    If argument left blank, loads whole data file
    '''
    def load_data(self,num=None,*cols):
        filename = self.data()
        data_cols = np.genfromtxt(filename,names=True,delimiter=',').dtype.names
        if cols == ():
            req_cols = data_cols
        else:
            req_cols = cols
        
        data_file_len = len(np.genfromtxt(filename,names=True,delimiter=','))
        
        if num == None:
            skip = 0
        elif num >= data_file_len:
            skip = 0
        else:
            skip = data_file_len - num
        
        bad_cols = []
        for col in cols:
            if col not in data_cols:
                bad_cols.append(col)
        if bad_cols != []:
            raise ValueError('{} not found for data column names.\n\nAvailable columns:\n{}'.format(bad_cols, data_cols))
            return_data = np.nan
        else:
            return_data = np.genfromtxt(filename,
                                        names=True,
                                        delimiter=',',
                                        usecols=req_cols,
                                        skip_footer=skip)
        return return_data
    
    '''
    returns path of any data that might have been written that day
    '''
    def data(self):
        if os.path.isfile(self.filename):
            return self.filename
        else:
            'There is no output data for {} yet.  Do write_data() first.'.format(self.data_path)


    '''
    playing around with different values of Rv
    '''
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
    
    #implementing the Deich Optical Photometry Extinction (DOPE) method
    #r1, r2 are the two colors F336W-F475W and F475W-F814W respectively
    def dope(self,r1,r2):

        intersect_function = self.findIntersection(self.IsochroneFit,
                                                   self.extinction_func,
                                                   r1,r2)
                                                   
        intersectpoint = fsolve(intersect_function,-1.5)
        newx = intersectpoint
        newy = self.IsochroneFit(intersectpoint)    
        dx = r1 - newx
        dy = r2 - newy

        a475 = dx/self.E336m475
        
        return (newx,newy,a475)
    
    
    #Creating the Panchromatic Hubble Andromeda Reddening Treasury (PHART)
    def phart(self,star):
        np.seterr(all='ignore')
        new336475,new475814,a475 = self.dope((star['F336W_VEGA']-star['F475W_VEGA']),
                                             (star['F475W_VEGA']-star['F814W_VEGA']))
        mag475 = star['F475W_VEGA']
        new475 = mag475 - 24.4 - a475


        return (new336475,new475814,new475,a475)
    

    def metal_fit(self,data_color,data_mag):
        z = self.metal_interp_function(data_color,data_mag)
        return (z)

    def write_data(self,lim=None):
        if lim is None:
            limit = len(self.raw_data)
        else:
            limit = lim
        
        import sys
        n = 0
        data_length = limit+0.0
        with open(self.filename, 'wb') as outfile:
            if self.is_sim == True:
                outfile.write('RA,DEC,F336W_NAKED-F475W_NAKED,F475W_NAKED-F814W_NAKED,F475W_NAKED,F336W-F475W_VEGA,F475W-F814W_VEGA,F475W_VEGA,F336W-F475W,F475W-F814W,F475W,A(F475W)_IN,A(F475W)_OUT,A(F475W)_DIFF,Z\n')
            elif self.is_sim == False:
                outfile.write('RA,DEC,F336W-F475W_VEGA,F475W-F814W_VEGA,F475W_VEGA,F336W-F475W,F475W-F814W,F475W,A(F475W),Z\n')
        with open(self.filename, 'ab') as outfile:
            for star in self.raw_data[:limit]:
                new336475, new475814, new475, a475 = self.phart(star)
                new336 = new336475 + new475
                new814 = (-1*new475814) - new475
                 
                z = self.metal_fit(new336475,new475)   
                if self.is_sim == True:
                    data = (star['RA'],
                            star['DEC'],
                            star['F336W_NAKED']-star['F475W_NAKED'],
                            star['F475W_NAKED']-star['F814W_NAKED'],
                            star['F475W_NAKED'],
                            star['F336W_VEGA']-star['F475W_VEGA'],
                            star['F475W_VEGA']-star['F814W_VEGA'],
                            star['F475W_VEGA'],
                            new336475[0],
                            new475814[0],
                            new475[0],
                            star['AF475W_IN'],
                            a475[0],
                            star['AF475W_IN']-a475[0],
                            z[0])        
                elif self.is_sim == False:
                    data = (star['RA'],
                            star['DEC'],
                            star['F336W_VEGA']-star['F475W_VEGA'],
                            star['F475W_VEGA']-star['F814W_VEGA'],
                            star['F475W_VEGA'],
                            new336475[0],
                            new475814[0],
                            new475[0],
                            a475[0],
                            z[0])
                outfile.write('{}\n'.format(','.join(map(str,(data)))))
                n += 1
                sys.stdout.write('\rDoing line {}, {}% done'.format(n,(n/data_length)*100))
                sys.stdout.flush()
      
class simulation(object):
    def __init__(self,numstars,iso_age):
        self.iso_age = iso_age
        self.iso_path = 'isochrones/{}Myr_finez.fits'.format(self.iso_age)
        self.isochrones = fits.open(self.iso_path)
        self.isochrones = self.isochrones[1].data
        self.isochrones = [self.isochrones[np.where(np.logical_and(self.isochrones["M_ini"]<10,
                                                                   self.isochrones["M_ini"]>2))]][0]
        self.iso_colors = np.rec.array([(self.isochrones["W336MAG"]-self.isochrones["F475MAG"]),
                                        (self.isochrones["F475MAG"]-self.isochrones["F814MAG"])],
                                        names=("Iso(336-475)","Iso(475-814)"))
                                        
        self.iso_color_x = sorted(self.iso_colors["Iso(336-475)"])
        self.iso_color_y = sorted(self.iso_colors["Iso(475-814)"])
        self.numstars = numstars
        self.mass_distribution = self.masslist(numstars)      
        
        self.month = time.strftime('%m')
        self.day = time.strftime('%d')
        self.filename = 'sim_out_{}_{}.fits'.format(self.month,self.day)
        
        self.write_data()       
 
    def masslist(self,numstars):
        population = np.random.random(numstars)
        intconst = 1.35
        normconst = 1.3527
        masses = (intconst/normconst)*(population**(1/(-intconst)))
        masses = masses[masses <= 50]
        return masses
    
    def interp(self,field):
        intp_function = interpolate.interp1d(self.isochrones['M_INI'],
                                      self.isochrones['{}'.format(field)],
                                      bounds_error = False)
        return intp_function
        
    def add_av(self):
        naked_w336 = self.interp('W336MAG')(self.mass_distribution)
        naked_f475 = self.interp('F475MAG')(self.mass_distribution)
        naked_f814 = self.interp('F814MAG')(self.mass_distribution)
        naked_w336f475 = naked_w336 - naked_f475
        naked_f475f814 = naked_f475 - naked_f814
        
        av_range = (3.5)*np.random.random(len(naked_f475))

        redvecs = {'w336f475':[0.446], 'f475f814':[0.61]}
        
        new_f475 = naked_f475-av_range
        dx = av_range*(redvecs['w336f475'])
        dy = av_range*(redvecs['f475f814'])
        new_w336f475 = naked_w336f475+dx
        new_f475f814 = naked_f475f814+dy
        
        new_w336 = new_w336f475 + new_f475
        new_f814 = -1 * (new_f475f814 - new_f475) 
        
        
        #ra, dec are dummy counters.  it's helpful for keeping track of individual stars
        #later and it makes it fit into the format that armpit expects.
        ra = np.linspace(1,len(naked_f475),len(naked_f475))
        dec = ra
        
        #armpit expects there to be an error column, so I use the following for all three.
        dummy = np.zeros(len(naked_f475))
        
        #I don't know how to prepare these columns for a FITS table file nicely,
        #so this is how it's going for now
        
        sim_data = fits.BinTableHDU.from_columns(
            [fits.Column(name='RA',format='E',array = ra),
             fits.Column(name='DEC',format = 'E',array = dec),
             fits.Column(name='MASS',format = 'E',array = self.mass_distribution),
             fits.Column(name='F336W_NAKED',format = 'E',array = naked_w336),
             fits.Column(name='F475W_NAKED',format = 'E',array = naked_f475),
             fits.Column(name='F814W_NAKED',format = 'E',array = naked_f814),
             fits.Column(name='F336W_VEGA',format = 'E',array = new_w336),
             fits.Column(name='F475W_VEGA',format = 'E',array = new_f475),
             fits.Column(name='F814W_VEGA',format = 'E',array = new_f814),
             fits.Column(name='AF475W_IN',format = 'E',array = av_range),
             fits.Column(name='F336W_ERR',format = 'E',array = dummy),
             fits.Column(name='F475W_ERR',format = 'E',array = dummy),
             fits.Column(name='F814W_ERR',format = 'E',array = dummy)])
        
        return sim_data
    
    def write_data(self):
        simulation_data = self.add_av()
        prihdr = fits.Header()
        prihdr['IS_SIM'] = True
        prihdr['ISO_AGE'] = self.iso_age
        prihdu = fits.PrimaryHDU(header=prihdr)
        
        hdulist = fits.HDUList([prihdu,simulation_data])
        
        hdulist.writeto('sim_out_{}_{}.fits'.format(self.month,self.day),clobber=True)