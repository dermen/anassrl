import os

from pylab import *
import h5py
import yaml

from popi.corr import correlate as Corr
from RingData import DiffCorr

def do_cors( p , ang, norm=True , trim = True):
    num_shot = p.shape[0]
    num_phi  = p.shape[1]
    if norm:
        for i in xrange( num_shot):
            p[i] /= ma.masked_equal( p[i], -1 ).mean()
        p [ p < 0] = -1
    
    c = zeros( ( num_shot-ang,num_phi))
    for i in xrange( num_shot - ang ):
        c[i] = Corr( p[i]-p[i+ang], p[i]-p[i+ang],0 )
    return c

ang = 16 # number of capillary rotations between shots that are subtracted
q = 10
polar_data = yaml.load(open('polar_imgs.yaml'))

cors = []
for run in polar_data:
    print run
    polar_imgs = np.load( '%s.npy'%polar_data[run]['polar_images'])[:,q].astype(float32)
    polar_mask = polar_data[run]['polar_mask'][q].astype(bool)
    polar_imgs[:,~polar_mask] = - 1
    
    correlations = do_cors( polar_imgs, ang, trim=False) 
    
    cors.append( correlations )

np.save( 'all_correlations', cors)

