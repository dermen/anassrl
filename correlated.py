import h5py
import yaml
from pylab import *

from RingData import DiffCorr

data_info     = yaml.load( open( 'silver_ssrl_data_info_final.yaml') )
data_hdf_indx = yaml.load( open( 'interped-shots_hdf5_index_mapping.yaml' ) )

data_hdf       = h5py.File( 'interped_shots.hdf5' )

cors =[]
for run in data_info:
    print run
    
    run_grp       = data_hdf[ run ]
    print run

    shot_pairs    = data_info[run]['good_exposure_pairs']
 
    path_to_indx  = ( (data_hdf_indx[run][s1], data_hdf_indx[s2])  for s1,s2 in shot_pairs  )
    indx          = [ (run_grp[p1], run_grp[p2]) for p1,p2 in path_to_indx]  

    data        = run_grp[ 'polardata'].value
    
    imshow( data[0], aspect='auto',vmax=10000)
    show()
    shot_diffs    = vstack( (  data_hdf[i1 ]-data[i2] for i1,i2 in indx ) )

    DC = DiffCorr( data )
    cor = DC.autocorr()
    cors.append( cor )

np.save( 'nearest_cors', cors ) 
