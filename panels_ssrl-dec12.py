from sys import argv
from scipy.interpolate import RectBivariateSpline as RBS
import numpy as np
import h5py
import matplotlib.cm as cm
from pylab import *
from popi.popi import polar
import os,re
#""" given a RBS `spline` object, evalute it at `interp_points`"""

# ~Some Functions
def spline_eval(spline, interp_points ):
    interped_vals = []
    for i in xrange( interp_points.shape[0] ):
        interped_vals.append(  spline( *interp_points[i] )[0][0] )
    interped_vals= array ( interped_vals )
    return interped_vals

def common_mode ( panel  ):
    """not used, but basic common_mode algorithm that doesnt work well"""
    intens_hist = np.histogram( panel.flatten(), bins = 10000 )
    intens_offset = (intens_hist[1][ argmax(intens_hist[0])] +\
                    intens_hist[1][ argmax(intens_hist[0]) +1  ] ) / 2.
    panel -= intens_offset
    return panel

def depolarize(panel,r_vals,phi_vals, wavelen,pixsize,detdist, pol_factor):
    """basic depolarization algorithm taken from
    J. Chem. Phys., Vol. 113, No. 20, 22 November 2000"""
    theta_vals = arctan2( r_vals*pixsize ,detdist  ) 
    norm  = pol_factor*( 1 - (sin(theta_vals)**2) * (cos( phi_vals )**2) )
    norm  += (1-pol_factor) * ( 1 - (sin(theta_vals)**2) * (sin( phi_vals )**2) )
    panel /= norm
    return panel

def find_nearest(array,value):
    """given an `array` of numbers, find the array index whose element
    is closest to `value`"""
    idx = (np.abs(array-value)).argmin()
    return idx

def window_rms(a1, window_size):
    a2 = power(a1,2)
    window = ones(window_size)/float(window_size)
    return sqrt(convolve(a2, window, 'valid'))

def get_rings(  bin_files, start, num_shots_per_proc, 
                qi,ai,bi, 
                panels, pan_inds, 
                outfilename, wavelen, pixsize, 
                detdist, pol_fac, num_phi ):
    """
    This is the work horse. Interpolates detector panels to form Bragg rings

    params
    ------
    bin_files, array-like
        list of bin file names
    
    start, int
        the starting bin file index

    num_shots_per_proc, int
        how many bin files to read after start
    
    qi, float
        initial guess at the radial position of ring (pixel units)
    
    ai,bi
        initial guess at the x,y beam center (pixel units)
        NOTE: x,y when looking at the bin image with deault  pylab.imshow
    
    panels, np.array
        multi dim array that defines the panel positions relative to one another
    
    pan_inds, np.array
        defines which panels are intersected by the Bragg ring
    
    outfileame, str
        name of the output file
    
    wavelen, float
        wavelength of beam (angstroms)
    
    pixsize, float
        size of pixel (assumed square pixels) (meters)
    
    detdist, float
        sample-detector distance (meters)
    
    pol_fac, float
        fraction of polarization in plane vs out of plane of beam traj
    
    num_phi, int
        number of interpolated pixels around the ring

    """
   
#   ~Form output arrays
    pol_intens     = ones( ( num_shots_per_proc, num_phi) )*-1
    pol_intens_bin = ones_like( pol_intens )*-1


    for i_shot in xrange( num_shots_per_proc):
#       ~Load image and find center
        
        print bin_files[ start + i_shot], i_shot

        img = fromfile( bin_files[ start + i_shot ],dtype=float32).reshape( (2463, 2527) )

        q = qi
        a = ai
        b = bi

#       ~Define full detector coordinates
        #X,Y        = meshgrid( arange(img.shape[0]), arange(img.shape[1]) ) # x,y dim of image
        X,Y        = meshgrid( arange(img.shape[1]), arange(img.shape[0]) ) # x,y dim of image
        
        R          = sqrt( (X-a)**2 + ( Y-b )**2 ) # distance of each pixel from center
        PHI        = arctan2( ( Y-b) , ( X-a ) )   # azimuthal angle of each pixel
        phi_values = linspace( -pi, pi, num_phi  ) # azimuthal angle of each interpolated pixel in each ring
#       ...and define x,y of each interpolated pixel around the ring
        ring       = array( [ [ q*sin(phi)+b, q*cos( phi) + a] for phi in phi_values ] )  
#       ~Do a panel-wise interpolation of the ring
        for j_pan, i_pan in pan_inds:
#           parameter "p" is [ymin,ymax,xmin,xmax] of panel            
            p    = panels[j_pan,i_pan]
            pan  = img[ p[0]:p[1], p[2]:p[3] ]
           # pan  = common_mode( pan)
            pan  = depolarize(  pan, R[ p[0]:p[1], p[2]:p[3] ], 
                                PHI[ p[0]:p[1], p[2]:p[3]    ], 
                                wavelen, pixsize,detdist, pol_fac )

#           ~Build the spline (rectBivariateSpline)
            rbs  = RBS( arange( p[0],p[1] ), arange( p[2],p[3]), pan,kx=5,ky=5 ) 
#           ...determine the x,y of each bragg-ring-pixel on this panel 
            interp_here =        ring[        ring[:,1] > p[2]   ]
            interp_here = interp_here[ interp_here[:,1] < p[3]-1 ]
            interp_here = interp_here[ interp_here[:,0] > p[0]   ]
            interp_here = interp_here[ interp_here[:,0] < p[1]-1 ]
            if list(interp_here) == []:
                print "moving on"
                continue
#           ...and compute the interpolated values
            
            vals = spline_eval( rbs, interp_here )
#           
            #~These are the phi_vals of the interpolated pixels
            phis = np.arctan2( interp_here[:,0]-b, interp_here[:,1]-a )
#           ...and the corresponding indices in our output array
            phi_inds = array( map( lambda x : find_nearest( phi_values, x ) , phis ) )
            pol_intens[ i_shot , phi_inds] = copy(vals)
            
            
            vm = vals.mean()
            vs = vals.std()
#           ~Make binary intensity
            #vals[vals > vm + 3.5*vs] = 0
            #vals[vals < vm - 3.5*vs] = 0

            vmask = ma.masked_equal( vals, 0 )


            #cutoff = vmask.mean() + vmask.std()
            cutoff = vm + vs
            vals[vals <= cutoff] = 0
            vals[ vals > 0]      = 1
            pol_intens_bin[ i_shot , phi_inds] = copy(vals)
           
            
#           ~Make some plots
            #subplot(212);
            #plot( phi_inds[ where( vals==1)[0] ], pol_intens[i_shot,phi_inds[where( vals==1 )[0]] ] ,'bd',lw=2)
            #plot( phi_inds, pol_intens[i_shot,phi_inds], 'bx',ms=2.25, lw=2)
        #xlabel('phi (0-2PI)',fontsize=22)
        #ylabel('counts',fontsize=22)
        #subplot(211)
        #imshow( img, vmin=0, vmax=5500 )
        #colorbar()
        #gca().add_patch( Circle(  (a,b), q, fc='none', lw=2, ec='k' ) )
        #show()

#   ~Output the rings
    out_file = h5py.File(outfilename, 'w' )
    out_file.create_dataset( 'pol_intens',     data=pol_intens )
    out_file.create_dataset( 'pol_intens_bin', data=pol_intens_bin )
    out_file.create_dataset( 'file_names', data=bin_files )
    out_file.close()


#=========================================================================================

# ~Expermimental parameters
wavelen = 0.7293   #angstrom
pixsize = 0.000172 #meter
#detdist = 0.275  #meter
#detdist = 0.245  #meter
pol_fac = 0.99

detdistances = { '188' : .188, '245' : .245,'275' : .275 }  #meter

detdist = detdistances[ argv[3] ]

# ~Pilatus 6m detector panels (modules) are on a 12 x 5 grid (60 panels total)
# ... we will provide the reference coordinates of each panel

# ~This is the y-coordinate of the top left corner of each panel when using default pylab.imshow( bin.data ) 
pan_ymin = array([0,    212,  424,  636, 
                  848,  1060, 1272, 1484,
                  1696, 1908, 2120, 2332])
# ~This is the x-coordinate of the top left corner of each panel when using default pylab.imshow( bin.data ) 
pan_xmin  = array([0, 494, 988, 1482, 1976])

# ~Dimensions of each panel (pixel unit)
pan_ydim  = 195
pan_xdim  = 487

# ~Mask this many pixels from each panel edge
mask_edg  = 3
pan_ymin += mask_edg
pan_xmin += mask_edg
pan_ydim = pan_ydim - mask_edg
pan_xdim = pan_xdim - mask_edg

# ~Store the 60 panel boundaries in [ymin,ymax,xmin,xmax] format
panels = zeros(  (pan_ymin.shape[0], pan_xmin.shape[0], 4 ) )
for j in xrange( pan_ymin.shape[0]):
    for i in xrange( pan_xmin.shape[0] ):
        panels[j,i,0] = pan_ymin[j]
        panels[j,i,1] = pan_ymin[j] + pan_ydim 
        panels[j,i,2] = pan_xmin[i]
        panels[j,i,3] = pan_xmin[i] + pan_xdim 

# ~Now provide a file with a list of bin file names
# ...each file should contain data on one exposure

#bin_files = map( lambda x:x.strip(), open( "binList.txt" ).readlines() )

# ~Beam center and Bragg ring parameters:

#qi      = 476.#925.5    #q of the Bragg ring in pixels
#ai      = 1231.#1231.85 #beam x_center in pixels
#bi      = 1264.5#1264.55 #beam y_center in pixels

#qi   = 591.1
#ai = 1231.0
#bi = 1264.45

#qi      = 925.5    #q of the Bragg ring in pixels
#ai      = 1231.85 #beam x_center in pixels
#bi      = 1264.55 #beam y_center in pixels

if argv[4] == '1':
# SILVER DET 188
    if argv[3] == '188':
        qi = 351.95
        ai = 1231.25
        bi = 1264.55
# SILVER DET 275
    elif argv[3] == '275':
        qi = 514.7
        ai = 1231.35
        bi = 1263.85

#SILVER DET 245
    elif argv[3] == '245':
        qi = 458.45
        ai = 1231.5
        bi = 1264.05

if argv[4] == '2':
# SILVER DET 188
    if argv[3] == '188':
        qi = 411.5
        ai = 1231.25
        bi = 1264.55
# SILVER DET 275
    elif argv[3] == '275':
        qi = 602.7
        ai = 1231.35
        bi = 1263.85

#SILVER DET 245
    elif argv[3] == '245':
        qi = 536.45
        ai = 1231.5
        bi = 1264.05

#num_phi = 2520 #6120   #num pixels to interpolate around Bragg ring

num_phi = 4320
# ~Now we must provide which panels
# the ring is intersecting (one could get fancy here and make this auto)

"""
panel_inds = [[5,1],
              [4,1],
              [3,2],
              [4,3],
              [5,3],
              [6,3],
              [7,3],
              [8,2],
              [7,1],
              [6,1]]

panel_inds = [[5,1],
              [4,1],
              [3,1],
              [3,2],
              [3,3],
              [4,3],
              [5,3],
              [6,3],
              [7,3],
              [8,3],
              [8,2],
              [8,1],
              [7,1],
              [6,1]]

"""
if argv[4] == '1':
#SILVER RING1 DET 188
    if argv[3] == '188':
        panel_inds = [[5,1],
                      [4,1],
                      [4,2],
                      [4,3],
                      [5,3],
                      [6,3],
                      [7,3],
                      [7,2],
                      [7,1],
                      [6,1]]

#SILVER RING1 DET 275
    if argv[3] == '275':
        panel_inds = [[5,1],
                      [4,1],              
                      [3,1],              
                      [3,2],              
                      [3,3],              
                      [4,3],              
                      [5,3],              
                      [6,3],              
                      [7,3],              
                      [8,3],              
                      [8,2],              
                      [8,1],              
                      [7,1],              
                      [6,1]]              

#SILVER RING1 DET 245
    if argv[3] == '245':
        panel_inds = [[5,1],
                      [4,1],
                      [4,2],
                      [3,2],
                      [4,3],
                      [5,3],
                      [6,3],
                      [7,3],
                      [7,2],
                      [8,2],
                      [7,1],
                      [6,1]]

if argv[4] == '2':
#SILVER RING1 DET 188
    if argv[3] == '188':
        panel_inds = [[5,1],
                      [4,1],
                      [4,2],
                      [4,3],
                      [5,3],
                      [6,3],
                      [7,3],
                      [7,2],
                      [7,1],
                      [6,1]]

#SILVER RING1 DET 275
    if argv[3] == '275':
        panel_inds = [[5,1],
                      [4,1],              
                      [3,1],              
                      [3,2],              
                      [3,3],              
                      [4,3],              
                      [5,3],              
                      [6,3],              
                      [7,3],              
                      [8,3],              
                      [8,2],              
                      [8,1],              
                      [7,1],              
                      [6,1]]              

#SILVER RING1 DET 245
    if argv[3] == '245':
        panel_inds = [[5,1],
                      [4,1],
                      [3,1],
                      [3,2],
                      [3,3],
                      [4,3],
                      [5,3],
                      [6,3],
                      [7,3],
                      [8,3],
                      [8,2],
                      [8,1],
                      [7,1],
                      [6,1]]
"""
panel_inds = [[ 5,0 ],
              [ 4,0 ],
              [ 3,0 ],
              [ 3,1 ],
              [ 2,1 ],
              [ 1,1 ],
              [ 1,2 ],
              [ 1,3 ],
              [ 2,3 ],
              [ 3,3 ],
              [ 3,4 ],
              [ 4,4 ],
              [ 5,4 ],
              [ 6,4 ],
              [ 7,4 ],
              [ 8,4 ],
              [ 8,3 ],
              [ 9,3 ],
              [10,3 ],
              [10,2 ],
              [10,1 ],
              [ 9,1 ],
              [ 8,1 ],
              [ 8,0 ],
              [ 7,0 ],
              [ 6,0 ]]
"""

# ~Apply a mask to the panels with detector shadow
# ..this is on every pilatus image we use from 12-2)
# ..and more specific masks can be applied to the interpolated pixels
for j,i in [[5,2],[5,3],[5,4] ]:
    panels[j,i,1] = panels[j,i,1] - 35
for j,i in [[6,2],[6,3],[6,4] ]:
    panels[j,i,0] = panels[j,i,0] + 40

if os.path.isdir( argv[1] ):
    bindir = argv[1]
    bin_files = [ x for x in  os.listdir(bindir ) if x.endswith('.bin')   ]
    bin_files = map( lambda x: os.path.join(bindir, x), bin_files ) 

elif os.path.isfile(argv[1]):
    bin_files = [os.path.abspath( argv[1] )]
# ~Vary these three params to parallelize the code
# ->file index in bin_files to start with 
start = 0
# ->number of files to continue analyzing after start
num_shots_per_proc = len(bin_files)
# ->name of the output file

outfilename = argv[2]

get_rings(  bin_files,start, num_shots_per_proc, qi,ai,bi, 
            panels,panel_inds, outfilename, wavelen, 
            pixsize, detdist, pol_fac, num_phi )
