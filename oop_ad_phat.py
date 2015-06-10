from astropy.io import fits
import numpy as np
import scipy.interpolate as interpolate
from scipy.optimize import fsolve


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
        
        if self.data_path.endswith('.csv'):
            self.raw_data = np.genfromtxt(self.data_path, names=True, delimiter=',')
        elif self.data_path.endswith('.fits'):
            self.raw_data = fits.open(self.data_path)
            self.raw_data = self.raw_data[1].data
            self.raw_data = self.raw_data[np.where(self.raw_data['F814W_ERR']<0.25)]
            self.raw_data = self.raw_data[np.where(self.raw_data['F475W_ERR']<0.25)]
        else:
            raise ValueError('Data must be in csv or fits formats')        
        

            
        #self.data_colors = np.rec.array([(self.raw_data['F336W_VEGA']-self.raw_data['F475W_VEGA']),
                                        #(self.raw_data['F475W_VEGA']-self.raw_data['F814W_VEGA'])],
                                        #names = ("Data(336-475)","Data(475-814)"))

        

        #self.raw_x = self.data_colors["Data(336-475)"]
        #self.raw_y = self.data_colors["Data(475-814)"]

        self.fit_coeffs = np.polyfit(self.iso_color_x,self.iso_color_y,3)

    '''
    Arguments: which columns you want, list
    Returns: data from inputted columns
    If argument left blank, loads whole data file
    '''
    def load_data(self,*cols):
        filename = self.data()
        data_cols = np.genfromtxt(filename,names=True,delimiter=',').dtype.names
        if cols == ():
            req_cols = data_cols
        else:
            req_cols = cols
            
        bad_cols = ()
        for col in cols:
            if col not in data_cols:
                bad_cols.append(col)
        if bad_cols != ():
            raise ValueError('{} not found for data names.  Available names:{}'.format(bad_cols, data_cols))
            return_data = np.nan
        else:
            return_data = np.genfromtxt(filename,names='True',delimiter=',',usecols=req_cols)
        return return_data
    
    '''
    returns path of any data that might have been written that day
    '''
    def data(self):
        import os.path
        import time
        month = time.strftime('%m')
        day = time.strftime('%d')
        filename = 'armpit_out/armpit_out_{}_{}.csv'.format(month,day)
        
        return filename if os.path.isfile(filename) else 'There is no data for {} yet.  Do write_data() first.'.format(self.data_path)


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
        initial_guesses = np.zeros(len(r1))-1.5
        
        intersect_function = self.findIntersection(self.IsochroneFit,
                                                   self.extinction_func,
                                                   r1,r2)
                                                   
        intersectpoint = fsolve(intersect_function,initial_guesses)
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

    

    def write_data(self):
        import time
        month = time.strftime('%m')
        day = time.strftime('%d')
        filename = 'armpit_out/armpit_out_{}_{}.csv'.format(month,day)
        
        import sys
        n = 0
        data_length = len(self.raw_data)+0.0
        with open(filename, 'wb') as outfile:
            outfile.write('RA,DEC,F336W-F475W_RAW,F475W-F814W_RAW,F475W_RAW,F336W-F475W,F475W-F814W,F475W,A(F475W),Z\n')
        
        
        with open(filename, 'ab') as outfile:
            for star in self.raw_data:
                new336475, new475814, new475, a475 = self.phart(star)
                new336 = new336475 + new475
                new814 = (-1*new475814) - new475
                
                z = self.metal_fit(new336475,new475)                
                data = (star['RA'],
                        star['DEC'],
                        star['F336W_VEGA']-star['F475W_VEGA'],
                        star['F475W_VEGA']-star['F814W_VEGA'],
                        star['F475W_VEGA'],
                        new336475,
                        new475814,
                        new475,
                        a475,
                        z)
                outfile.write('{}\n'.format(','.join(map(str,(data)))))
                n += 1
                sys.stdout.write('\rDoing line {}, {}% done'.format(n,(n/data_length)*100))
                sys.stdout.flush()
      
class simulation(object):
    def __init__(self,iso_age,numstars):
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
 
    def masslist(self,numstars):
        population = np.random.random(numstars)
        intconst = 1.35
        normconst = 1.3527
        masses = (intconst/normconst)*(population**(1/(-intconst)))
        masses = masses[masses <= 50]
        return masses
    
    def interp(self,field):
        intp_function = intp.interp1d(self.isochrones['M_INI'],
                                      self.isochrones['{}'.format(field)],
                                      bounds_error = False)
        return intp_function
        
    def add_grid(self):
        naked_w336 = self.interp('W336MAG')(self.mass_distribution)
        naked_f475 = self.interp('F475MAG')(self.mass_distribution)
        naked_f814 = self.interp('F814MAG')(self.mass_distribution)
        naked_w336f475 = naked_w336 - naked_f475
        naked_f475f814 = naked_f475 - naked_f814
        
        av_range = np.linspace(0,3.5,num=200)  