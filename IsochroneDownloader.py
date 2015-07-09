import numpy as np
import scipy.interpolate as intp
import mechanize
import re
import os



# Downloads Padova isochrones of many ages
# with the parameters I want and converts them to fits
# mechanize is difficult to do correctly. It is almost definitely not
# totally kosher in this implementation.

isochrone_dir = '.' #INSERT ISOCHRONE DIRECTORY HERE

for age in np.arange(5,101,.5):
    response = mechanize.urlopen('http://stev.oapd.inaf.it/cgi-bin/cmd')
    forms = mechanize.ParseResponse(response)
    response.close()
    iso_form = forms[0]
    filter_menu = iso_form.find_control('photsys_file')
    for item in filter_menu.items:
        if item.name == 'tab_mag_odfnew/tab_mag_wfc3_wide.dat':
            item.selected = True

    iso_form['isoc_val'] = ['2']
    iso_form['isoc_age0'] = '{}e6'.format(age)
    request2 = iso_form.click()
    response2 = mechanize.urlopen(request2)
    output_page = response2.read()
    p = re.compile('output\d+\.dat')
    iso_filename = p.search(output_page).group(0)
    dat_filename = '{}Myr_finez.dat'.format(age)
    fits_filename = '{}Myr_finez.fits'.format(age)    
    br = mechanize.Browser()
    br.retrieve(
                'http://stev.oapd.inaf.it/~lgirardi/tmp/{}'.format(iso_filename),
                os.path.join('./{}'.format(isochrone_dir),dat_filename.format(age))
                )
    
    
    os.system('./stilts tcopy {}/{} {}/{} ifmt=ascii ofmt=fits'.format(isochrone_dir,
                                                                       dat_filename,
                                                                       isochrone_dir,
                                                                       fits_filename))
    os.system('rm {}/{}'.format(isochrone_dir,
                                dat_filename))

    print '{}Myr isochrone downloaded & converted'.format(age)
        
