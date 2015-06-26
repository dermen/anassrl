import h5py
import os
from pylab import *
import numpy as np

from RingData import DiffCorr

#files1 = ['188-%d.hdf5'%x for x in xrange( 1,19)]
#files1 += ['275-%d.hdf5'%x for x in xrange( 1,11)]
#files1 += ['245-%d.hdf5'%x for x in xrange( 0,4)]

#files1 = map( lambda x:  os.path.abspath(x), files1 ) 
#files1 = filter( lambda x: os.path.exists(x) , files1)
#hdfs1 = [h5py.File( f,'r' ) for f in files1]

dat_dir = '/Users/mender/data_analysis/analysis_output'
os.chdir(dat_dir)

files1 = [ x for x in os.listdir('.') if x.endswith('polar.npy') ]
pols1 = [ load( pf).astype(float32) for f in files1 ]
num_pols = len(pols1)

from popi.corr import correlate as Corr

def do_cors( p1 ,p2, ang, norm=True ):
    num_shot = p1.shape[0]
    assert( p1.shape[1] == p2.shape[1] )
    num_phi = p1.shape[1]
    if norm:
        for i in xrange( num_shot):
            p1[i] /= ma.masked_equal( p1[i], -1 ).mean()
            p2[i] /= ma.masked_equal( p2[i], -1 ).mean()
        p1 [ p1 < 0] = -1
        p2 [ p2 < 0] = -1

    c = zeros( ( num_shot-ang,num_phi))
    for i in xrange( num_shot - ang ):
        c[i] = Corr( p1[i]-p1[i+ang], p2[i]-p2[i+ang],0 )
    return c

ang = 16
cor = [  do_cors( pols1[x][ order[x] ], pols1[x][ order[x] ], ang ) for x in xrange( num_pols )  ]


out = h5py.File( 'diff_cors_%d_fnames-auto.hdf5'%ang, 'w')

for i in xrange( num_pols ) :
    out.create_dataset( 'cors/cor%d'%i, data=cor[i] )

out.create_dataset( 'scans', data=array( files) )
out.close()
