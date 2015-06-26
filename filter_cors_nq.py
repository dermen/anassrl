from sys import argv
from scipy import optimize
from pylab import *
import numpy as np
import h5py

def find_nearest(array,value):
    """given an `array` of numbers, find the array index whose element
    is closest to `value`"""
    idx = (np.abs(array-value)).argmin()
    return idx

#q_ind = 10

cor = np.load(  argv[1] )
cor = vstack( [ c.astype(float32) for c in cor ] )

num_sig = 3.5 # SSRL
num_bin = cor.shape[0] / 10.
edge =  200


def curve_fitfun( X,p0,p1, p2  ):
    return (p2/pi) * (p1 / 2. ) * (1/( (X-p0)*(X-p0) + p1*p1/4. )   )

total_del = 0

num_phi = cor.shape[1]

for ang in xrange( edge, num_phi / 2 - edge ) :
    H = np.histogram( cor[:,ang],bins=num_bin,normed=True )
    y = H[0]
    x = [ ( H[1][i] + H[1][i+1] )/ 2. for i in xrange( H[1].shape[0]-1) ]

    amp = y.max()
    mu = x[argmax(y)]
    sig = amp/1000.
    parm = [ mu+1e-6, sig, amp]
    fit,success = optimize.curve_fit( curve_fitfun,xdata=x,ydata=y,p0=tuple(parm) )
    
    print mu, sig, amp
    mu  = fit[0] 
    sig = fit[1]
    amp = fit[ 2 ] 
    
    cutoff_H = mu + num_sig * sig*3
    cutoff_L = mu - num_sig * sig*3
    
#    plot( x, y )
#    plot ( x, curve_fitfun( x, *fit), lw=2, c='r' ) 
#    plot( ones(2)* cutoff_H, [-1,amp],'--',c='lime',lw=3,alpha=0.6)
#    plot( ones(2)* cutoff_L, [-1,amp],'--',c='lime',lw=3,alpha=0.6)
#    draw()

    bad = where( logical_or( cor[:,ang] > cutoff_H, cor[:,ang] < cutoff_L ))[0]
   
    good = list ( set(range(cor.shape[0] ) ) - set(bad) )

#    print >> fnames_out, " ".join( fnames[ bad ] ) 

    #cor = delete( cor, bad,axis = 0 )

    cor = cor[good]
#   fnames = fnames[ good ]

    #plot( x, fitfun(fit,x), 'r',lw=2)
    #ylim(0,y.max() )
    
    #hist( cor[:,ang],bins=num_bin,normed=True,log=True)
    #show()

    total_del += len(bad)
    print "total deleted shots:",total_del, 'remaining angles:', num_phi/2 - edge - ang

#fnames_out.close()

np.save( argv[1] + 'filtered', cor )



