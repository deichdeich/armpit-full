import reddening as reddening

print "\n##########\n\nGrabbing extinction\n\n############\n\n"
reddening.mp_write()



import metal_fitter as mf
print "\n#########\n\nInterpolating for metallicity\n\n############\n\n"
mf.write_to_file()

execfile('do_plots.py')
# be able to pass data, isochrones as arguments
