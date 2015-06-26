import os

import json
import h5py
from pylab import *

from RingData import RingFit, InterpSimple


def norm_polar_images( imgs, mask_val = -1 ): 
    norms = ma.masked_equal( imgs, mask_val).mean(axis=2)
    imgs /= norms[:,:,None]
    imgs[ imgs < 0 ] = mask_val
    return imgs

data_info = json.load( open( 'silver_ssrl_data_info_final.json', 'r') )


# DICTIONARY DESCRIPTION (data_info)
# ==================================
# Each key of this dictionary is a string
# unique to a set of exposures (a run).

# FIn this experiment, a
# run corresponds to a 
# rotation of the capillary
# about the azimuth, as the x-ray beam 
# exposing. 

# About every 3 exposures ,the exposed sample
# regions become uncorrelated

# Each key unlocks a sub-dictionary with the following keys:
# ========================================================
# 'exposure_time_seconds',
# 'type',
# 'wavelength_angstroms',
# 'pixel_size_meters',
# 'detector_distance_meters',
# 'filenames_on_texas',
# 'file_names_on_drive_mac'
# 'slow_fast_dimensions_Cstyle',
# 'good_exposure_pairs'

centering_q = 2.66 # where to look in recirpocal space
num_phi     = 4320 # pixels around each inteprolated ring
mask_file   = 'ssrl_mask_5.bin' # bool; same structure as images, 0 is masked, 1 is kept
nq          = 3 # extent of polar image from  q_111 in pixel units
prefix      = 'interped_shots'

output_hdf = h5py.File( '%s.hdf5'%prefix, 'w' )

fname_to_hdf_map = {}

for run in data_info:
    print ('run: %s'%run)
    run_info   = data_info[run]
#   load the run meta data
    img_shape  = run_info['slow_fast_dimensions_Cstyle']
    wavelen    = run_info['wavelength_angstroms' ]
    pixsize    = run_info['pixel_size_meters']
    detdist    = run_info['detector_distance_meters']
    filenames  = run_info[ 'filenames_on_texas' ] 
    datatype   = run_info['type']
 
    num_imgs   = len(filenames)
    
#   some useful functions
    pix2invang  = lambda qpix : sin( arctan( qpix*pixsize/detdist ) /2 )* 4 * pi / wavelen 
    invang2pix  = lambda qia  : tan(2*arcsin(qia*wavelen/4/pi))*detdist/pixsize
    
#   center the first image (remove this if you already know the center)
    test_img      = fromfile( filenames[0], dtype=datatype ).reshape( img_shape)
    y_init        = img_shape[0]/2. # slow dim
    x_init        = img_shape[1]/2. # fast dim (last axis in numpy array)
    q_init        = invang2pix(centering_q)
    
    RF = RingFit( test_img )
    x_center, y_center, q_ring = RF.fit_circle( beta_i=(x_init, y_init,q_init) )
    
#   initialize the interpolator
    interpolater = InterpSimple( x_center, y_center, 
                           q_ring+nq, q_ring-nq,
                           num_phi,raw_img_shape=img_shape)
    I_ = interpolater.nearest

#   load a mask file( same shape as the images, but boolean)
    mask_img    = fromfile( mask_file, dtype=bool).reshape(img_shape)
    pmask       = I_( mask_img, dtype=bool).round()
    output_hdf.create_dataset( '%s/polar_mask'%run ,data = pmask.astype(int))

#   generates the images 
    img_gen    = ( fromfile(f, dtype=datatype).reshape(img_shape) for f in filenames)

#   iterate over the images : interpolate and save them to disk
    polar_imgs  = array( [ pmask* I_( img_gen.next()) for dummie in xrange(num_imgs)] )

    polar_imgs = norm_polar_images( polar_imgs , mask_val = 0) 

    output_hdf.create_dataset( '%s/polar_data'%run ,data=  polar_imgs)

#   save a lookup-map
    for indx,f in enumerate( filenames):
        fname_to_hdf_map[f] = ( run, indx )

#   save meta data
    output_hdf.create_dataset( '%s/x_center'%run,   data = x_center)
    output_hdf.create_dataset( '%s/y_center'%run,   data = y_center)
    output_hdf.create_dataset( '%s/q_ring'%run,     data = q_ring)
    output_hdf.create_dataset( '%s/num_phi'%run,     data = num_phi)
    output_hdf.create_dataset( '%s/wavelen'%run ,   data = wavelen)
    output_hdf.create_dataset( '%s/pixsize'%run ,   data = pixsize)
    output_hdf.create_dataset( '%s/detdist'%run ,   data = detdist)

#   save the q-mapping
    q_in_pix = arange( q_ring-nq, q_ring + nq )
    q_map    = array( [ [ ind, pix2invang(q) ] for ind,q in enumerate( q_in_pix) ] )
    output_hdf.create_dataset( '%s/q_mapping'%run ,    data = q_map )

# save the data
output_hdf.close()
# save the map
dump_f = open( '%s_hdf5_indx_map.json'%prefix, 'w')
json.dump( fname_to_hdf_map, dump_f )
dump_f.close()


