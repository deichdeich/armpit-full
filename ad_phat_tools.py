'''
ad_phat_tools.py:
alex deich PHAT analysis tools.

Author: Alex Deich
Date: June 12 2015
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
import sys

class armpit(object):

    def __init__(self,data_path,iso_age,extinct_col = None,
    usefile = None,rv = 3.1,filenum = 0,comments = ''):
    
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
        self.extinct_col = extinct_col

        # get the path of the fits file 
        ### TO DO: make it accept csv files   
        if type(self.data_path) is StellarPopulationSimulator:
            self.raw_file = self.data_path.output_file
        elif type(self.data_path) is str:
            if self.data_path.endswith('.fits'):
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
            self.sim_params = self.data_header['FREEPRMS']
        else:
            self.is_sim = False
            self.sim_params = np.nan
        
        # trim the data
        self.raw_data = self.raw_fits[1].data
        self.raw_data = self.raw_data[np.where(self.raw_data['F814W_ERR']<0.25)]
        self.raw_data = self.raw_data[np.where(self.raw_data['F475W_ERR']<0.25)]


        # calculate the polynomial fit of the isochrone
        self.fit_coeffs = np.polyfit(self.iso_color_x,self.iso_color_y,3)


        self.month = time.strftime('%m')
        self.day = time.strftime('%d')
        self.usefile = usefile
        self.comments = comments
        self.Rv = rv
        self.E336m475,self.E475m814 = self.redvec(self.Rv)
        if self.usefile == None:
            self.filename = 'multi_age_analysis/armpit_out_{}_{}_{}_{}.csv'.format(self.iso_age,
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
        data_cols = np.genfromtxt(filename,names=True,delimiter=',',skip_header = 4).dtype.names
        if cols == ():
            req_cols = data_cols
        else:
            req_cols = cols
        
        data_file_len = len(np.genfromtxt(filename,names=True,delimiter=',',skip_header = 4))

        
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
                                        skip_header = 4,
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
    def redvec(self,Rv_value):
        if Rv_value == 3.1:
            return (0.446,0.610)
        elif Rv_value == 5:
                return (0.208,0.467)
        else:
                return ValueError("{} is not a valid Rv value".format(Rv_value))
    


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
        if new336475 < 2.0 and new475814<-0.6:
            new336475,new475814,a475 = (np.nan,np.nan,np.nan)
        mag475 = star['F475W_VEGA']
        
        # 24.4: andromeda distance modulus
        intrinsic_475 = mag475 - 24.4 - a475

        if ~np.isnan(new336475) or ~np.isnan(new475814) or ~np.isnan(intrinsic_475) or ~np.isnan(a475):
            return (new336475[0],new475814[0],intrinsic_475[0],a475[0])
        else:
            return (np.nan,np.nan,np.nan,np.nan)
    
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
        csvcomments = '# Simulation: {}\n'.format(self.is_sim)+'# Free parameters in simulation: {}\n'.format(self.sim_params)+'# Fitted Age: {}\n'.format(self.iso_age)+'# Comments: {}\n'.format(self.comments)
        if self.is_sim == True:
            csvhdr = csvcomments+'RA,DEC,MASS,F336W_NAKED-F475W_NAKED,F475W_NAKED-F814W_NAKED,F475W_NAKED,F336W-F475W_VEGA,F475W-F814W_VEGA,F475W_VEGA,F336W-F475W,F475W-F814W,F475W,AV_IN,AV,AV_DIFF,Z\n'
            
        elif self.is_sim == False:
            csvhdr = csvcomments+'RA,DEC,F336W-F475W_VEGA,F475W-F814W_VEGA,F475W_VEGA,F336W-F475W,F475W-F814W,F475W,AV,Z\n'
        return csvhdr
        
       
    def write_data(self,lim=None,draw_plots=False):
        if lim is None:
            limit = len(self.raw_data)
        else:
            limit = lim
        
     
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
                            new336475,
                            new475814,
                            new475,
                            star['AV_IN'],
                            a475,
                            star['AV_IN']-a475,
                            z)        
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
                sys.stdout.write('\rDoing line {}, {}% done     '.format(n,(n/data_length)*100))
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
        if point_list == {}:
            self.point_list = point_list
        else:
            self.point_list = pickle.load(open(point_list,'rb'))
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

class simulation(object):

    # Integrating dN/dM ~ M^(-2.35) to get masses of a population
    def masslist(self,numstars):
        population = np.random.random(numstars)
        intconst = 1.35
        normconst = 1.3527
        masses = (intconst/normconst)*(population**(1/(-intconst)))
        return masses
# to do: make it grab the two closest isochrones, not just rounding to
# the nearest integer, so that you can have any range of isochrones in
# the directory.
class ArtificialInterpolatedMagnitudes(simulation):
    def __init__(self):
        self.interp_dict = {}
        self.iso_dir = 'shitload_of_isochrones/'
        self.filters = ['F336W','F475W','F814W']
        self.iso_dict = {}
        
        print '\nCollecting isochrones...'
        self.fill_iso_dict()
        print 'Found {} isochrones with age range from {}Myr to {}Myr.'.format(len(self.iso_dict),min(self.iso_dict.keys()),max(self.iso_dict.keys()))
        
    def __call__(self,age,mass,z_val,solar_norm = False):
        self.mass = mass
        self.age = float(age)
        self.z_val = z_val
        self.age_gap = self.get_age_gap()
        self.solar_norm = solar_norm
        self.iso_names = {'upper':'{}'.format(self.age_gap[0]),
                          'lower':'{}'.format(self.age_gap[1])}
        
        mags = self.simulate()
        return mags
        #print 'Magnitude in F336W: {}\nMagnitude in F475W: {}\nMagnitude in F814W: {}'.format(mags[0],mags[1],mags[2])
    
    
    # looks at available isochrones and fills iso_dict with isochrones and their
    # corresponding ages
    def fill_iso_dict(self):
        iso_list = os.listdir(self.iso_dir)
        self.fits_files = [self.iso_dir+file for file in iso_list if '.fits' in file]
        if self.fits_files ==[]:
            raise IOError('Isochrones must be in fits format')
        else:
            for file in self.fits_files:
                isochrone = fits.open(file)[1].data
                if 'log(age/yr)' in isochrone.dtype.names:
                    iso_age = (10**float(isochrone['log(age/yr)'][1]))/1e6
                    self.iso_dict[iso_age] = file
                else:
                    raise ValueError('Isochrones must have \'log(age/yr)\' column as provided by Padova.  Offending file: {}'.format(file))
                 
    
    # looks at iso_dict and finds the closest two isochrones
    def get_age_gap(self):
        available_ages = np.zeros(len(self.fits_files))
        counter=0
        for key in self.iso_dict:
            available_ages[counter] = key
            counter+=1
        
        available_ages = np.sort([available_ages])[0]
        nearest_isochrone_index = (np.abs(available_ages-self.age)).argmin()
        nearest_isochrone_age = available_ages[nearest_isochrone_index]
        if nearest_isochrone_age > self.age:
            if nearest_isochrone_age == np.min(available_ages):
                self.upper_age = self.lower_age = nearest_isochrone_age
                raise RuntimeWarning('{} younger than maximum available isochrone age of {}.  This may lead to an unreliable simulated value.'.format(self.age, nearest_isochrone_age))
                return (self.iso_dict[nearest_isochrone_age],self.iso_dict[nearest_isochrone_age])
            else:
                next_isochrone_age = available_ages[nearest_isochrone_index-1]
                self.upper_age = nearest_isochrone_age
                self.lower_age = next_isochrone_age
                return (self.iso_dict[nearest_isochrone_age],self.iso_dict[next_isochrone_age])
        elif nearest_isochrone_age < self.age:
            if nearest_isochrone_age == np.max(available_ages):
                self.upper_age = self.lower_age = nearest_isochrone_age
                raise RuntimeWarning('{} older than maximum available isochrone age of {}.  This may lead to an unreliable simulated value.'.format(self.age, nearest_isochrone_age))
                return (self.iso_dict[nearest_isochrone_age],self.iso_dict[nearest_isochrone_age])
            else:
                next_isochrone_age = available_ages[nearest_isochrone_index+1]
                self.upper_age = next_isochrone_age
                self.lower_age = nearest_isochrone_age
                return (self.iso_dict[next_isochrone_age],self.iso_dict[nearest_isochrone_age])
        
        
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

        lower_weight = (self.age-self.lower_age)/(self.upper_age-self.lower_age)
        upper_weight = 1-lower_weight

        f336w_sim = (upper_weight*uf336w_sim)+(lower_weight*lf336w_sim)
        f475w_sim = (upper_weight*uf475w_sim)+(lower_weight*lf475w_sim)
        f814w_sim = (upper_weight*uf814w_sim)+(lower_weight*lf814w_sim)
        
        return f336w_sim,f475w_sim,f814w_sim
    
    def interpolate(self,iso_file):
        isochrone = fits.open(iso_file)[1].data
        mass = isochrone['M_ini'].copy()
        if self.solar_norm == False:
            metals = isochrone['Z'].copy()
        elif self.solar_norm == True:
            metals = isochrone['Z'].copy()
            metals = np.log10(metals/0.019)
        else:
            raise ValueError('solar_norm must be True or False, not {}.'.format(self.solar_norm))

        points = zip(mass,metals)
        filt_dict = {}
        for filt in self.filters:
            filt_dict[filt] = interpolate.LinearNDInterpolator(points,isochrone[filt].copy())
        
        self.interp_dict[iso_file] = filt_dict

class StellarPopulationSimulator(simulation):
    def __init__(self,numstars,ages='all',z_values='all',masses='all',av='none',rv=3.1,output_file='',fits_comments=''):
        self.numstars = int(numstars)
        self.z_values = z_values
        self.masses = masses
        self.av = av
        self.rv = rv
        self.ages = ages
        self.month = time.strftime('%m')
        self.day = time.strftime('%d')
        self.fits_comments = fits_comments
        self.output_file = output_file
        if self.output_file == '':
            self.output_file = 'sim_pops/sim_out_{}_{}_{}.fits'.format('%.2E'%self.numstars,
                                                                self.month,
                                                                self.day)
        else:
            if type(self.output_file) is not str:
                raise TypeError('Any user-defined output filename must be a string')
        
        self.simstar = ArtificialInterpolatedMagnitudes()
        
        self.age_lim = (max(self.simstar.iso_dict.keys()),
                        min(self.simstar.iso_dict.keys()))
        
        
        
        self.free_parameters = []
        
        if numstars == '':
            if len(ages) == len(z_values) == len(masses):
                numstars = len(ages)
        
            else:
                raise ValueError('All input parameters must have the same length')
        
        if self.ages == 'all':
            self.free_parameters.append('ages')
            #self.age_list = (100-4)*np.random.random(self.numstars)-4
        else:
            try:
                if len(ages) == self.numstars:
                    self.age_list = ages
                elif len(ages) == 1:
                    self.age_list = np.zeros(self.numstars)+age[0]
                elif len(ages) != self.numstars or len(self.ages) != 1:
                    raise ValueError('Age array must equal number of stars or 1.')
                elif max(ages) > self.age_lim[0] or min(ages) < self.age_lim[1]:
                    raise ValueError('Maximum allowed age: {}Myr, minimum allowed age: {}Myr'.format(age_lim[0],age_lim[1]))
            except TypeError:
                self.age_list = np.zeros(self.numstars)+ages
        
        if self.z_values == 'all':
            self.free_parameters.append('z')
            #self.z_list = (0.03-0.0001)*np.random.random(self.numstars)-0.0001
        else:
            try:
                if len(z_values) == self.numstars:
                    self.z_list = z_values
                elif len(z_values) == 1:
                    self.z_list = np.zeros(self.numstars)+self.z_values[0]
                elif len(z_values) != self.numstars or len(self.z_values) != 1:
                    raise ValueError('Z array must equal number of stars or 1.')
                elif max(z_values) > 0.03 or min(self.z_values) < 0.0001:
                    raise ValueError('Maximum allowed Z: 0.03, minimum allowed Z: 0.0001')
            except TypeError:
                self.z_list = np.zeros(numstars)+z_values
        
        if self.masses == 'all':
            self.free_parameters.append('masses')
            #self.mass_list = masslist(numstars)
        else:
            try:
                if len(masses) == numstars:
                    self.mass_list = masses
                elif len(masses) == 1:
                    self.mass_list = np.zeros(numstars)+mass[0]
                elif len(masses) != numstars or len(masses) != 1:
                    raise ValueError('Mass array must equal number of stars or 1.')
            except TypeError:
                self.mass_list = np.zeros(numstars)+masses
        
        if self.av == 'none':
            self.av_list = np.zeros(self.numstars)
        elif self.av == 'grid':
            self.av_list = (3.5)*np.random.random(self.numstars)
        else:
            try:
                if len(av) == numstars:
                    self.av_list = av
                elif len(av) == 1:
                    self.av_list = np.zeros(numstars)+av[0]
                elif len(av) != numstars or len(av) != 1:
                    raise ValueError('Av array must equal number of stars or 1.')
            except TypeError:
                self.av_list = np.zeros(numstars)+av
                
        if self.av is not 'none' and self.rv == 3.1:
            self.ext_vals = {'f336w':1.667,'f475w':1.221,'f814w':0.61}
        elif self.av is not 'none' and self.rv == 5:
            self.ext_vals = {'f336w':1.349,'f475w':1.14,'f814w':0.67}
        elif self.av is not 'none':
            raise ValueError('Rv must be either 3.1 or 5.  Not {}.'.format(self.rv))
        
        #initialize all the arrays
        self.ra_arr = np.zeros(self.numstars)
        self.dec_arr = np.zeros(self.numstars)
        self.mass_arr = np.zeros(self.numstars)
        self.age_arr = np.zeros(self.numstars)
        self.z_arr = np.zeros(self.numstars)
        self.f336w_naked_arr = np.zeros(self.numstars)
        self.f475w_naked_arr = np.zeros(self.numstars)
        self.f814w_naked_arr = np.zeros(self.numstars)
        self.av_arr = np.zeros(self.numstars)
        self.f336w_arr = np.zeros(self.numstars)
        self.f475w_arr = np.zeros(self.numstars)
        self.f814w_arr = np.zeros(self.numstars)
        self.dummy_arr = np.zeros(self.numstars)
        self.starlist()
        print '\ndone with simulation, writing FITS...'
        self.write_to_file()
        print '\nFITS file {} written'.format(self.output_file)
    
    def simulate(self,starnum):
        perc_done = round(25*(starnum/float(self.numstars)))
        progress_bar = '=' * int(perc_done+1)
        sys.stdout.write('\rSimulating star {}... |{}{}| {}%'.format(starnum,
                                                            progress_bar,
                                                             ' '*(25-int(perc_done+1)),
                                                             round(100*((starnum+1)/float(self.numstars)),2)))
        sys.stdout.flush()
        self.age,self.mass,self.z,self.av = self.get_params(starnum)
        self.star = [0,0,0]
    
        if self.free_parameters != []:
            while (self.star[0]-self.star[1]>-0.5 or np.isnan(self.star[0])
                    or np.isnan(self.star[1]) or np.isnan(self.star[2])): 
                self.age,self.mass,self.z,self.av = self.get_params(starnum)
                self.star = self.simstar(self.age,self.mass,self.z)
        elif self.free_parameters == []:
            self.star = self.simstar(self.age,self.mass,self.z)
        
        self.mass_arr[starnum] = self.mass
        self.age_arr[starnum] = self.age
        
        self.ra_arr[starnum] = self.dec_arr[starnum] = starnum
        
        self.f336w = self.star[0]
        self.f475w = self.star[1]
        self.f814w = self.star[2]
    
        self.f336w_naked_arr[starnum] = self.f336w
        self.f475w_naked_arr[starnum] = self.f475w
        self.f814w_naked_arr[starnum] = self.f814w
        
        self.a336 = self.av*self.ext_vals['f336w']
        self.a475 = self.av*self.ext_vals['f475w']
        self.a814 = self.av*self.ext_vals['f814w']
    
        self.av_arr[starnum] = self.av

        # add av
        self.f336w += self.a336
        self.f475w += self.a475
        self.f814w += self.a814
            
        self.f336w_arr[starnum] = self.f336w
        self.f475w_arr[starnum] = self.f475w
        self.f814w_arr[starnum] = self.f814w
        self.z_arr[starnum] = self.z
    
    
    
    def get_params(self,starnum):    
        if 'ages' in self.free_parameters:
            self.star_age = (self.age_lim[0]-self.age_lim[1])*np.random.random(1)+self.age_lim[1]
        elif 'ages' not in self.free_parameters:
            self.star_age = self.ages[starnum]
        
        if 'masses' in self.free_parameters:
            self.star_mass = self.masslist(1)
        elif 'masses' not in self.free_parameters:
            self.star_mass = self.mass_list[starnum]
            
        if 'z' in self.free_parameters:
            self.z_val = (0.03-0.0001)*np.random.random(1)+0.0001
        elif 'z' not in self.free_parameters:
            self.z_val = self.z_list[starnum] 
        
        self.star_av = self.av_list[starnum]  
    
        return self.star_age,self.star_mass,self.z_val,self.star_av
        
    def starlist(self):
        stars = xrange(self.numstars)
        map(self.simulate,stars)
        
    def get_fits_data(self):
        sim_data = fits.BinTableHDU.from_columns(
            [fits.Column(name='RA',format='E',array = self.ra_arr),
             fits.Column(name='DEC',format = 'E',array = self.dec_arr),
             fits.Column(name='AGE',format='E',array = self.age_arr),
             fits.Column(name='MASS',format = 'E',array = self.mass_arr),
             fits.Column(name='Z',format = 'E',array = self.z_arr),
             fits.Column(name='F336W_NAKED',format = 'E',array = self.f336w_naked_arr),
             fits.Column(name='F475W_NAKED',format = 'E',array = self.f475w_naked_arr),
             fits.Column(name='F814W_NAKED',format = 'E',array = self.f814w_naked_arr),
             fits.Column(name='F336W_VEGA',format = 'E',array = self.f336w_arr),
             fits.Column(name='F475W_VEGA',format = 'E',array = self.f475w_arr),
             fits.Column(name='F814W_VEGA',format = 'E',array = self.f814w_arr),
             fits.Column(name='AV_IN',format = 'E',array = self.av_arr),
             fits.Column(name='F336W_ERR',format = 'E',array = self.dummy_arr),
             fits.Column(name='F475W_ERR',format = 'E',array = self.dummy_arr),
             fits.Column(name='F814W_ERR',format = 'E',array = self.dummy_arr)])
        return sim_data
    
    def write_to_file(self):
        sim_data = self.get_fits_data()
        prihdr = fits.Header()
        prihdr['IS_SIM'] = True
        prihdr['FREEPRMS'] = ','.join(self.free_parameters)
        if self.av == 'none':
            self.rv = 'no extinction added'
        prihdr['RV_VALUE'] = self.rv
        prihdr['COMMENTS'] = self.fits_comments
        prihdu = fits.PrimaryHDU(header=prihdr)
        hdulist = fits.HDUList([prihdu,sim_data])
        hdulist.writeto(self.output_file,clobber=True)