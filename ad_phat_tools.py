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
import os
import time
import matplotlib.pyplot as plt
import numpy.lib.recfunctions as rfn
import pickle

class armpit(object):

    def __init__(self,data_path,iso_age,extinct_col = None,
    usefile = None,filenum = 0):
    
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
        self.filenum = filenum
        # get the path of the fits file 
        
        ### TO DO: make it accept csv files   
        if type(self.data_path) is SingleAgeSimulation:
            self.raw_file = self.data_path.output_file
        elif self.data_path.endswith('.fits'):
            self.raw_file = self.data_path
        else:
            raise IOError('Data input must either be a fits file containing the relevant columns or a simulation object.')
        
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
        self.usefile = usefile
        if self.usefile == None:
            self.filename = 'multi_age_analysis/armpit_out_beastav_{}_{}_{}_{}.csv'.format(self.iso_age,
                                                                                        self.month,
                                                                                        self.day,
                                                                                        self.filenum)
        else:
            self.filename = self.usefile
        
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
    
    #If you want to de-extinct the star with some other av included in the input file
    def subtract_av(self,star,avcol):
        av = star[avcol]
        intrinsic_475 = star['F475W_VEGA'] - 24.4 - av
        dx = av*self.E336m475
        dy = av*self.E475m814
        
        new336475 = (star['F336W_VEGA']-star['F475W_VEGA'])-dx
        new475814 = (star['F475W_VEGA']-star['F814W_VEGA'])-dy

        return(new336475,new475814,intrinsic_475,av)        

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
                if self.extinct_col == None:
                    new336475, new475814, new475, a475 = self.phart(star)
                else:
                    new336475, new475814, new475, a475 = self.subtract_av(star,self.extinct_col)
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
                            new336475,
                            new475814,
                            new475,
                            a475,
                            z)
                outfile.write('{}\n'.format(','.join(map(str,(data)))))
                n += 1
                sys.stdout.write('\rDoing line {}, {}% done'.format(n,(n/data_length)*100))
                sys.stdout.flush()
        if draw_plots == True:
            self.make_plot('cmd')
            print 'doing cmd...'
            self.make_plot('ccd')
            print 'doing ccd...'
            
class region_draw(object):
### To do: be able to write out and save a dictionary of regions
###  and then be able to pass it this dictionary again.
###  and figure out why the cmd is wonky
###  and make all plots look final

    def __init__(self,
                    ax,
                    map_fig, 
                    zplot_fig,
                    cmdplot_fig,
                    ccdplot_fig,
                    avplot_fig,
                    cmd_ax,
                    ccd_ax,
                    avhist_ax,
                    zhist_ax,
                    data,
                    point_list = {}):
        self.previous_point = []
        self.start_point = []
        self.end_point = []
        self.line = None    
        self.point_list = point_list
        self.map_fig = map_fig
        self.map_fig.canvas.draw()
        self.map_fig.canvas.set_window_title('Spatial Map')
        self.data = data
        self.zplot_fig = zplot_fig
        self.zplot_fig.canvas.set_window_title('Metallicity Histogram')
        self.avplot_fig = avplot_fig
        self.avplot_fig.canvas.set_window_title('Av Histogram')
        self.cmdplot_fig = cmdplot_fig
        self.cmdplot_fig.canvas.set_window_title('Color-Magnitude Diagram')
        self.ccdplot_fig = ccdplot_fig
        self.ccdplot_fig.canvas.set_window_title('Color-Color Diagram')
        self.cmd_ax = cmd_ax
        self.ccd_ax = ccd_ax
        self.avhist_ax = avhist_ax
        self.zhist_ax = zhist_ax
        self.ax = ax
        self.star_arr = np.zeros((len(self.data),1),dtype=[('REGION_NUM',int)])
        self.star_arr['REGION_NUM'] += 999
        self.region_arr = rfn.merge_arrays([self.data,self.star_arr],flatten = True)
        self.region_arr = np.ma.masked_array(self.region_arr,np.isnan(self.region_arr['Z']))
        
        self.region_counter = 0
        
        self.draw_spatial()
        
        self.color_list = ['red','blue','green','purple','black']
        self.color_num = 0
        self.color_mod = 0
        
        self.region_data = 0
        self.region_check = False
    def __call__(self):
        self.point_list = {}
    def button_press_callback(self, event):
        if event.inaxes:
        #if there is no previous_point, create a new key in the point_list dictionary.
        #in either case, put x,y in the most recent key.
            if self.previous_point == []:
                self.region_counter += 1
                self.point_list[self.region_counter] = []
            x, y = event.xdata, event.ydata
            if event.button == 1:  # If you press the right button
                    if self.line == None: # if there is no line, create a line
                        self.line = plt.Line2D([x,  x],
                                           [y, y],
                                           marker = 'o',
                                           color = self.color_list[self.color_num],
                                           linewidth = 2)
                        self.start_point = [x,y]
                        self.previous_point =  self.start_point 
                        self.ax.add_line(self.line)
                        self.map_fig.canvas.draw()
                    # add a segment
                    else: # if there is a line, create a segment
                        self.line = plt.Line2D([self.previous_point[0], x], 
                                           [self.previous_point[1], y],
                                           marker = 'o',
                                           color = self.color_list[self.color_num],
                                           linewidth = 2)
                        self.previous_point = [x,y]
                        event.inaxes.add_line(self.line)
                        self.map_fig.canvas.draw()
                
                    self.point_list[self.region_counter].append((x,y))
                    
            elif event.button == 3 and self.line != None: # close the loop
                        print self.point_list.keys()
                        self.line = plt.Line2D([self.previous_point[0],
                                            self.start_point[0]],
                                           [self.previous_point[1],
                                            self.start_point[1]],
                                            marker = 'o',
                                            color = self.color_list[self.color_num],
                                            linewidth = 2)
                        self.color_mod += 1
                        self.color_num = self.color_mod%len(self.color_list)
                        self.previous_point = []
                        self.ax.add_line(self.line)
                        self.map_fig.canvas.draw()
                        self.line = None
                        self.region_check = False
                        
                        
            
            
    def key_press_callback(self,event): 
        # the keys I want to define are: clear point_list and plots: (c)
        #                                make cmd,ccd,zhist,avhist: ('1','2','3','4')
        #                                save cmd,ccd,zhist,zvhist: (ctrl + 1,2,3 or 4)
        #                                do/save all: (a)/(ctrl+a) 
        
        
        if event.key == '1':
            self.cmdplot_fig.canvas.draw()
            if self.point_list == {}:
                color_data = self.data['F336WF475W_VEGA']
                mag_data = self.data['F475W_VEGA']-24.4
                z_data = np.log10(self.data['Z']/0.019)
            else:
                self.set_data()
                color_data = self.region_data['F336WF475W_VEGA']
                mag_data = self.region_data['F475W_VEGA']-24.4
                z_data = np.log10(self.region_data['Z']/0.019)
            
            self.make_cmd(self.cmd_ax,
                          self.cmdplot_fig,
                          color_data,
                          mag_data,
                          z_data)
        
        elif event.key == '2':
            self.ccdplot_fig.canvas.draw()
            if self.point_list == {}:
                color1_data = self.data['F336WF475W_VEGA']
                color2_data = self.data['F475WF814W_VEGA']
                av_data = self.data['AF475W'] 
            else:
                self.set_data()
                color1_data = self.region_data['F336WF475W_VEGA']
                color2_data = self.region_data['F475WF814W_VEGA']
                av_data = self.region_data['AF475W']
                    
            self.make_ccd(self.ccd_ax,
                          self.ccdplot_fig,
                          color1_data,
                          color2_data,
                          av_data)
        
        elif event.key == '3':
            self.zplot_fig.canvas.draw()
            if self.point_list == {}:
                histdata = self.data['Z']
                self.make_hist(self.zhist_ax,self.zplot_fig,histdata,1,-1,.25)
            else:
                self.set_data()
                for region in self.point_list:
                    histdata = self.region_data[self.region_data['REGION_NUM']==region]
                    histdata = np.log10(histdata['Z']/0.019)
                    self.make_hist(self.zhist_ax,self.zplot_fig,histdata,region,-1,.25)
                
                self.zhist_ax.set_title('Histogram of log(Z/0.019) for selected regions')
        
        elif event.key == '4':
            self.avplot_fig.canvas.draw()
            if self.point_list == {}:
                histdata = self.data['AF475W']
                self.make_hist(self.avhist_ax,self.avplot_fig,histdata,1,0,max(self.data['AF475W']))
            else:
                self.set_data()
                for region in self.point_list:
                    histdata = self.region_data[self.region_data['REGION_NUM']==region]
                    histdata = histdata['AF475W']
                    self.make_hist(self.avhist_ax,self.avplot_fig,histdata,region,0,max(self.data['AF475W']))
                
                self.avhist_ax.set_title('Histogram of A(F475W) for selected regions')
        
        elif event.key == 'ctrl+1':
            extent = self.cmd_ax.get_window_extent().transformed(self.cmdplot_fig.dpi_scale_trans.inverted())
            self.plot_fig.savefig('awesomecmd.png', bbox_inches = extent)
            print 'CMD saved '
        
        elif event.key == 'ctrl+2':
            extent = self.ccd_ax.get_window_extent().transformed(self.ccdplot_fig.dpi_scale_trans.inverted())
            self.plot_fig.savefig('awesomeccd.png', bbox_inches = extent)
            print 'CCD saved '
            
        elif event.key == 'ctrl+3':
            extent = self.zhist_ax.get_window_extent().transformed(self.zplot_fig.dpi_scale_trans.inverted())
            self.plot_fig.savefig('awesomezhist.png', bbox_inches = extent)
            print 'Metallicity histogram saved '   
        
        elif event.key == 'ctrl+4':
            extent = self.avhist_ax.get_window_extent().transformed(self.avplot_fig.dpi_scale_trans.inverted())
            self.plot_fig.savefig('awesomeavhist.png', bbox_inches = extent)
            print 'Av Histogram saved '
        
        elif event.key == 'w':
            ##write out point list
            month = time.strftime('%m')
            day = time.strftime('%d')
            pickle.dump(self.point_list,open('skyselect_regions_{}_{}.p'.format(month,day),'wb'))
        
        elif event.key == 'c':
            self.region_counter = 0
            self.point_list = {}
            self.cmd_ax.clear()
            self.ccd_ax.clear()
            self.zhist_ax.clear()
            self.avhist_ax.clear()
            self.ax.clear()
            self.draw_spatial()
            self.map_fig.canvas.draw()
            self.zplot_fig.canvas.draw()
            self.avplot_fig.canvas.draw()
            self.ccdplot_fig.canvas.draw()
            self.cmdplot_fig.canvas.draw()
            self.color_num = 0
            self.color_mod = 0
            self.line = None
            self.region_check = False
            self.region_arr['REGION_NUM'] = 999
            
    def set_data(self):
        if self.region_check == False:
            print 'Getting stars in region(s)\n'
            self.region_data = self.get_star_in_region()
            self.region_check = True
            for region in self.point_list:
                print 'found {} stars in region {}'.format(len(self.region_data[self.region_data['REGION_NUM'] == region]),region)
        else:
            pass
    
    def get_star_in_region(self):
        for star in self.region_arr:
            if star['REGION_NUM'] == 999:
                for region in self.point_list:
                    if self.point_in_poly(star['RA'],star['DEC'],self.point_list[region]):
                        star['REGION_NUM'] = region
            else:
                pass
        return self.region_arr[self.region_arr['REGION_NUM']!=999]
    
    def draw_spatial(self):
        self.ax.set_title('Spatial Map')
        self.ax.scatter(self.data['RA'],
                        self.data['DEC'],
                        s = 2,
                        edgecolors = 'none',
                        c=self.data['Z'])
    
    def make_cmd(self,ax,fig,magdata,colordata,z_data):
        ax.scatter(self.data['F336WF475W_VEGA'],self.data['F475W_VEGA']-24.4,edgecolors='none',c='gray')
        ax = ax.scatter(magdata,colordata,edgecolors='none',c=z_data)
        ax.set_title('CMD of selected region(s)')
        ax.set_xlim(-2,0)
        ax.set_ylim(-6,0)
        ax.invert_yaxis()
        fig.canvas.draw()
    def make_ccd(self,ax,fig,color1,color2,av_data):
        ax.scatter(self.data['F336WF475W_VEGA'],self.data['F475WF814W_VEGA'],edgecolors='none',c='gray')
        ax.scatter(color1,color2,edgecolors='none',c=av_data)
        ax.set_title('CCD of selected region(s)')
        ax.set_xlim(-2,0)
        ax.set_ylim(-1,2)
        cb = fig.colorbar()
        cb.set_clim(0,6)
        fig.canvas.draw()
    
    def make_hist(self,ax,fig,value,region_color,xmin,xmax):
        plot_data = value[~np.isnan(value)]
        histcolor = self.color_list[region_color-1]
        n, bins, patches = ax.hist(plot_data, 
                                   50, 
                                   normed=1,
                                   cumulative = True,
                                   histtype='step',
                                   color = histcolor)
        ax.set_xlim(xmin,xmax)
        fig.canvas.draw()
    
    def point_in_poly(self,x,y,poly):
        n = len(poly)
        inside = False
        p1x,p1y = poly[0]
        for i in xrange(n+1):
            p2x,p2y = poly[i % n]
            if y > min(p1y,p2y):
                if y <= max(p1y,p2y):
                    if x <= max(p1x,p2x):
                        if p1y != p2y:
                            xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or x <= xints:
                            inside = not inside
            p1x,p1y = p2x,p2y
        return inside

class armplot(object):
    def __init__(self,data_object):
        print '\n loading {}... \n'.format(data_object.filename)
    
        self.data = data_object.load_data()
        
        print '\n{} with {} data points loaded.'.format(data_object.filename,
                                                      len(self.data))
    def skyselect(self):
        map_fig = plt.figure()
        ax = map_fig.add_subplot(111)
        
        zplot_fig = plt.figure()
        cmdplot_fig = plt.figure()
        ccdplot_fig = plt.figure()
        avplot_fig = plt.figure()
        cmd_ax = cmdplot_fig.add_subplot(111)
        ccd_ax = ccdplot_fig.add_subplot(111)
        avhist_ax = avplot_fig.add_subplot(111)
        zhist_ax = zplot_fig.add_subplot(111)
        ax.set_title('')
        cursor = region_draw(ax,map_fig,zplot_fig,cmdplot_fig,ccdplot_fig,avplot_fig,cmd_ax,
                                                ccd_ax,
                                                avhist_ax,
                                                zhist_ax,self.data)
        map_fig.canvas.mpl_connect('button_press_event', cursor.button_press_callback)
        map_fig.canvas.mpl_connect('key_release_event',cursor.key_press_callback)
        zplot_fig.canvas.mpl_connect('key_release_event',cursor.key_press_callback)
        avplot_fig.canvas.mpl_connect('key_release_event',cursor.key_press_callback)
        ccdplot_fig.canvas.mpl_connect('key_release_event',cursor.key_press_callback)
        cmdplot_fig.canvas.mpl_connect('key_release_event',cursor.key_press_callback)
        plt.show()

class SingleAgeSimulation(object):
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
        
        new_f475 = naked_f475+av_range
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

class ArtificialInterpolatedMagnitudes(object):
    def __init__(self):
        self.interp_dict = {}
        self.iso_dir = 'shitload_of_isochrones'
        self.filters = ['F336W','F475W','F814W']
    
    def __call__(self,age,mass,z_val,av,solar_norm = False):
        self.mass = mass
        self.age = float(age)
        self.z_val = z_val
        self.av = av
        self.age_gap = self.get_age_gap()
        self.solar_norm = solar_norm
        self.iso_names = {'upper':'{}/{}Myr_finez.fits'.format(self.iso_dir,self.age_gap[0]),
                          'lower':'{}/{}Myr_finez.fits'.format(self.iso_dir,self.age_gap[1])}
       
        self.a336 = self.av*1.667
        self.a475 = self.av*1.219
        self.a814 = self.av*0.61
        
        
        mags = self.simulate()
        print 'Magnitude in F336W: {}\nMagnitude in F475W: {}\nMagnitude in F814W: {}'.format(mags[0],mags[1],mags[2])
    
    def get_age_gap(self):
        rounded_age = round(self.age)
        if rounded_age > self.age:
            return (rounded_age-1,rounded_age)
        elif rounded_age < self.age:
            return (rounded_age,rounded_age+1)
        else:
            return (self.age,self.age)
        
    def simulate(self):
        if self.iso_names['upper'] not in self.interp_dict:
            self.interpolate(self.iso_names['upper'])
        elif self.iso_names['upper'] in self.interp_dict:
            pass
        
        if self.iso_names['lower'] not in self.interp_dict:
            self.interpolate(self.iso_names['lower'])
        elif self.iso_names['lower'] in self.interp_dict:
            pass
        
        uf336w_sim,uf475w_sim,uf814w_sim = [self.interp_dict[self.iso_names['upper']][filt](self.mass,self.z_val) for filt in self.filters]
        lf336w_sim,lf475w_sim,lf814w_sim = [self.interp_dict[self.iso_names['lower']][filt](self.mass,self.z_val) for filt in self.filters]
        
        upper_weight = abs(self.age_gap[1]-self.age)
        lower_weight = abs(self.age_gap[0]-self.age)
        if upper_weight == lower_weight == 0:
            upper_weight = lower_weight = 0.5
        
        
        f336w_sim = (upper_weight*uf336w_sim)+(lower_weight*lf336w_sim)+self.a336
        f475w_sim = (upper_weight*uf475w_sim)+(lower_weight*lf475w_sim)+self.a475
        f814w_sim = (upper_weight*uf814w_sim)+(lower_weight*lf814w_sim)+self.a814
        
        return f336w_sim,f475w_sim,f814w_sim
    
    def interpolate(self,iso_file):
        isochrone = fits.open(iso_file)[1].data
        mass = isochrone['M_ini'].copy()
        if self.solar_norm == False:
            metals = isochrone['Z'].copy()
        elif self.solar_norm == True:
            metals = isochrone['Z'].copy()
            metals = np.log10(metals/0.019)

        points = zip(mass,metals)
        filt_dict = {}
        for filt in self.filters:
            filt_dict[filt] = interp.LinearNDInterpolator(points,isochrone[filt].copy())
        
        self.interp_dict[iso_file] = filt_dict
