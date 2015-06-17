'''
ad_phat_tools.py:
alex deich PHAT analysis tools.

Author: Alex Deich
Date: Jun 12 2015
'''

from astropy.io import fits
import numpy as np
import scipy.interpolate as interpolate
from scipy.optimize import fsolve
import os.path
import time
import matplotlib.pyplot as plt

class armpit(object):

    def __init__(self,data_path,iso_age):
    
        # get isochrone and chop out the part we're interested in (the nondegenerate part)
        self.iso_age = iso_age
        self.iso_path = 'isochrones/{}Myr_finez.fits'.format(self.iso_age)
        self.isochrones = fits.open(self.iso_path)
        self.isochrones = self.isochrones[1].data
        self.isochrones = [self.isochrones[np.where(np.logical_and(self.isochrones["M_ini"]<10,
                                                                   self.isochrones["M_ini"]>2))]][0]
                                                                   
        self.isochrones = [self.isochrones[np.where((self.isochrones['W336MAG']-self.isochrones['F475MAG'])<0)]][0]
                                                                   
        self.iso_colors = np.rec.array([(self.isochrones["W336MAG"]-self.isochrones["F475MAG"]),
            (self.isochrones["F475MAG"]-self.isochrones["F814MAG"])],
            names=("Iso(336-475)","Iso(475-814)"))
                                        
        self.iso_color_x = sorted(self.iso_colors["Iso(336-475)"])
        self.iso_color_y = sorted(self.iso_colors["Iso(475-814)"])
        
        self.iso_mag = self.isochrones['F475MAG']
        self.zvalues = self.isochrones['Z']
        
        #make sure you *don't* use the sorted color for the metallicity interpolation
        self.interp_points = zip(self.iso_colors['Iso(336-475)'],self.iso_mag)
        self.metal_interp_function = interpolate.LinearNDInterpolator(self.interp_points,
                                                                      self.zvalues)
        
        
        self.data_path = data_path
        
        # get the path of the fits file 
        
        ### TO DO: make it accept csv files   
        if type(self.data_path) is simulation:
            self.raw_file = self.data_path.output_file
        elif self.data_path.endswith('.fits'):
            self.raw_file = self.data_path
        else:
            raise IOError('Data input must either be a fits or csv file with containing the relevant columns or a simulation object.')
        
        self.raw_fits = fits.open(self.raw_file)
        
        self.data_header = self.raw_fits[0].header
        
        # check if it's a simulation output-- I *could* just do this when I do
        # if type(...) above, but I can conceive of having simulation output
        # that's not just an object.
        if 'IS_SIM' in self.data_header:
            self.is_sim = self.data_header['IS_SIM']
            self.sim_age = self.data_header['SIM_AGE']
        else:
            self.is_sim = False
            self.sim_age = np.nan
        
        # trim the data
        self.raw_data = self.raw_fits[1].data
        self.raw_data = self.raw_data[np.where(self.raw_data['F814W_ERR']<0.25)]
        self.raw_data = self.raw_data[np.where(self.raw_data['F475W_ERR']<0.25)]


        # calculate the polynomial fit of the isochrone
        self.fit_coeffs = np.polyfit(self.iso_color_x,self.iso_color_y,3)


        self.month = time.strftime('%m')
        self.day = time.strftime('%d')
        self.filename = 'armpit_out/armpit_out_{}myr_{}_{}.csv'.format(self.iso_age,
                                                                       self.month,
                                                                       self.day)
        
    '''
    Arguments: which columns you want, list
    Returns: data from inputted columns
    If argument left blank, loads whole data file
    '''
    def load_data(self,*cols):
        filename = self.data()
        data_cols = np.genfromtxt(filename,names=True,delimiter=',',skip_header=3).dtype.names
        if cols == ():
            req_cols = data_cols
        else:
            req_cols = cols
        
        data_file_len = len(np.genfromtxt(filename,names=True,delimiter=',',skip_header=3))

        
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
                                        skip_header=3,
                                        usecols=req_cols)
        return return_data
    
    '''
    returns path of any data that might have been written that day
    '''
    def data(self):
        if os.path.isfile(self.filename):
            return self.filename
        else:
            raise IOError('There is no output data for {} yet.  Do write_data() first.'.format(self.data_path))


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
        
        # 24.4: andromeda distance modulus
        intrinsic_475 = mag475 - 24.4 - a475


        return (new336475,new475814,intrinsic_475,a475)
    

    def metal_fit(self,data_color,data_mag):
        z = self.metal_interp_function(data_color,data_mag)
        return (z)


    def csv_header(self):
        csvcomments = '# Simulation: {}\n'.format(self.is_sim)+'# Simulated population age: {}\n'.format(self.sim_age)+'# Fitted Age: {}\n'.format(self.iso_age)
        if self.is_sim == True:
            csvhdr = csvcomments+'RA,DEC,MASS,F336W_NAKED-F475W_NAKED,F475W_NAKED-F814W_NAKED,F475W_NAKED,F336W-F475W_VEGA,F475W-F814W_VEGA,F475W_VEGA,F336W-F475W,F475W-F814W,F475W,A(F475W)_IN,A(F475W),A(F475W)_DIFF,Z\n'
            
        elif self.is_sim == False:
            csvhdr = csvcomments+'RA,DEC,F336W-F475W_VEGA,F475W-F814W_VEGA,F475W_VEGA,F336W-F475W,F475W-F814W,F475W,A(F475W),Z\n'
        return csvhdr
        
       
    def write_data(self,lim=None,draw_plots=False):
        if lim is None:
            limit = len(self.raw_data)
        else:
            limit = lim
        
        import sys
        n = 0
        data_length = limit+0.0
        with open(self.filename, 'wb') as outfile:
            outfile.write(self.csv_header())
            
        with open(self.filename, 'ab') as outfile:
            for star in self.raw_data[:limit]:
                new336475, new475814, new475, a475 = self.phart(star)
                new336 = new336475 + new475
                new814 = (-1*new475814) - new475
                 
                z = self.metal_fit(new336475,new475)   
                if self.is_sim == True:
                    data = (star['RA'],
                            star['DEC'],
                            star['MASS'],
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
        if draw_plots == True:
            self.make_plot('cmd')
            print 'doing cmd...'
            self.make_plot('ccd')
            print 'doing ccd...'
            
    # plotting method
    def make_plot(self,type):
        import matplotlib.pyplot as plt
        plottables = self.load_data('F336WF475W_VEGA',
                                    'F475WF814W_VEGA',
                                    'F475W_VEGA',
                                    'F336WF475W',
                                    'F475WF814W',
                                    'F475W',
                                    'AF475W',
                                    'Z')
        z_data = np.ma.masked_array(plottables,np.isnan(plottables['Z']))
        if type == 'ccd':
            plt.figure()
            plt.scatter(plottables['F336WF475W_VEGA'],
                        plottables['F475WF814W_VEGA'],
                        edgecolors='none',
                        color='gray')
            plt.scatter(z_data['F336WF475W_VEGA'],
                        z_data['F475WF814W_VEGA'],
                        edgecolors = 'none',
                        c=z_data['AF475W'])
            cb = plt.colorbar()
            
            cb.set_clim([-0.5,3.5])
            cb.set_label('A(F475W)')
            
            plt.scatter(self.iso_color_x,self.iso_color_y,color='black')
            
            plt.xlim(-2,0)
            plt.xlabel('F336W-F475W')
            plt.ylim(-1,2)
            plt.ylabel('F475W-F814W')
            plt.savefig('ccd{}.png'.format(self.iso_age))
            
        elif type == 'cmd':
            plt.figure()
            plt.scatter(plottables['F336WF475W_VEGA'],
                        plottables['F475W_VEGA'],
                        edgecolors='none',
                        color='gray')
            plt.scatter(z_data['F336WF475W_VEGA'],
                        z_data['F475W_VEGA']-24.4,
                        edgecolors = 'none',
                        c=np.log10(z_data['Z']/0.019))
            cb = plt.colorbar()

            
            #cb.set_clim([-0.5,3.5])
            cb.set_label('log(Z/0.019)')
            
            #plt.scatter(self.iso_color_x,self.iso_mag,color='black')
            
            plt.xlim(-2,0)
            plt.ylim(-6,0)
            plt.xlabel('F336W-F475W')
            plt.ylabel('F475W')
            ax = plt.gca()
            ax.invert_yaxis()
            plt.savefig('cmd{}.png'.format(self.iso_age))

# region_draw class adapted from code by Daniel Kornhauser
class region_draw(object):
    def __init__(self, ax, fig):
        self.previous_point = []
        self.start_point = []
        self.end_point = []
        self.line = None    
        self.point_list = []
        self.fig =  fig
        self.fig.canvas.draw()
        
    def button_press_callback(self, event):
        if event.inaxes: 
            x, y = event.xdata, event.ydata
            ax = event.inaxes
            if event.button == 1:  # If you press the right button
                    if self.line == None: # if there is no line, create a line
                        self.line = plt.Line2D([x,  x],
                                           [y, y],
                                           marker = 'x')
                        self.start_point = [x,y]
                        self.previous_point =  self.start_point 
                        ax.add_line(self.line)
                        self.fig.canvas.draw()
                    # add a segment
                    else: # if there is a line, create a segment
                        self.line = plt.Line2D([self.previous_point[0], x], 
                                           [self.previous_point[1], y],
                                           marker = 'x')
                        self.previous_point = [x,y]
                        event.inaxes.add_line(self.line)
                        self.fig.canvas.draw()
                    
                    self.point_list.append((x,y))
                    
            elif event.button == 3 and self.line != None: # close the loop
                        print self.point_list
                        self.line = plt.Line2D([self.previous_point[0],
                                            self.start_point[0]],
                                           [self.previous_point[1],
                                            self.start_point[1]],
                                            marker = 'x')
                 
                        ax.add_line(self.line)
                        self.fig.canvas.draw()
                        self.line = None
                        
        


class armplot(object):
    def __init__(self,):
        self.blah = 'blah'
    def skyselect(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('')
        cursor = region_draw(ax,fig)
        fig.canvas.mpl_connect('button_press_event', cursor.button_press_callback)
        plt.show()
    #def make_cmd(self,data):
    
    #def make_ccd(self,data):
    #def make_hist(self,data,col)
    

class simulation(object):
    def __init__(self,numstars,sim_age):
        self.sim_age = sim_age
        self.iso_path = 'isochrones/{}Myr_finez.fits'.format(self.sim_age)
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
        self.output_file = 'sim_pops/sim_out_{}myr_{}_{}.fits'.format(self.sim_age,
                                                                   self.month,
                                                                   self.day)
        
        self.write_data()       
 
    # Integrating dN/dM ~ M^(-2.35) to get masses of a population
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
    
    # adds a grid of Av
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
        prihdr['SIM_AGE'] = self.sim_age
        prihdu = fits.PrimaryHDU(header=prihdr)
        
        hdulist = fits.HDUList([prihdu,simulation_data])
        
        hdulist.writeto(self.output_file,clobber=True)
