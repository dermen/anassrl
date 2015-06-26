from scipy import optimize
from pylab import *
import numpy as np
import h5py

def find_nearest(array,value):
    """given an `array` of numbers, find the array index whose element
    is closest to `value`"""
    idx = (np.abs(array-value)).argmin()
    return idx

# filter data
from sys import argv

f = h5py.File(argv[1])
cor = vstack(( f['cors/'+'cor%d'%i].value for i in xrange( 32))    )

fnames = [ f['cors/'+'fnames%d'%i].value for i in xrange( 32 )    ]
fnames = vstack( fnames)
fnames = array( map( lambda x: x[0]+'---'+x[1], fnames  )) 


num_sig = 7.5 # SSRL
#num_sig = 2.5 # SACLA
num_bin = 1000
edge =  300

# gauss
#fitfun = lambda p,X:  p[0]*exp(-((X-p[1])**2)/(2*p[2]*p[2]))
# lorentz
#fitfun = lambda p,X:  p[0] * (1/pi) * (p[2] / 2. ) * (1/( (X-p[1])*(X-p[1]) + p[2]*p[2]/4. )   )
fitfun = lambda p,X:  (1/pi) * (p[1] / 2. ) * (1/( (X-p[0])*(X-p[0]) + p[1]*p[1]/4. )   )
errfun = lambda p,X,Y: fitfun(p,X) - Y

def curve_fitfun( X,p0,p1  ):
    return (1/pi) * (p1 / 2. ) * (1/( (X-p0)*(X-p0) + p1*p1/4. )   )


total_del = 0

num_phi = cor.shape[1]

for ang in xrange( edge, num_phi / 2 - edge ) :
    # filter using cor
    H = np.histogram( cor[:,ang],bins=num_bin,normed=True )
    y = H[0]
    x = [ ( H[1][i] + H[1][i+1] )/ 2. for i in xrange( H[1].shape[0]-1) ]

    #amp = y.max()
    mu = x[argmax(y)]
    sig = abs(x[find_nearest(y, 50)] - mu )
    parm = [ mu+1e-6, sig+1e-6]    
    fit,success = optimize.curve_fit( curve_fitfun,xdata=x,ydata=y,p0=tuple(parm) )

    mu = fit[0] 
    sig = fit[1]
    
    bad1 = where( logical_or( cor[:,ang] > mu + num_sig*sig, cor[:,ang] < mu - num_sig*sig ))[0]
   

#   filter using cor2
    #H = np.histogram( cor2[:,ang],bins=num_bin,normed=True )
    #y = H[0]
    #x = [ ( H[1][i] + H[1][i+1] )/ 2. for i in xrange( H[1].shape[0]-1) ]

    #amp = y.max()
    #mu = x[argmax(y)]
    #sig = abs(x[find_nearest(y, 50)] - mu )
    #parm = [ mu+1e-6, sig+1e-6]    
    #fit,success = optimize.curve_fit( curve_fitfun,xdata=x,ydata=y,p0=tuple(parm) )

    #mu = fit[0] 
    #sig = fit[1]
    
    #bad2 = where( logical_or( cor[:,ang] > mu + num_sig*sig, cor[:,ang] < mu - num_sig*sig ))[0]
   
    #bad = unique( array( list(bad1) + list(bad2) ) ) 

    good = list ( set(range(cor.shape[0] ) ) - set(bad1) )
    
    cor = cor[good]
    #cor2 = cor2[good]
    #cor3 = cor3[good]

    #fnames3 = fnames3[good]

    #plot( x, fitfun(fit,x), 'r',lw=2)
    #ylim(0,y.max() )
    
    #plot( ones(2)* (mu + 2*sig), [-1e3,1e3],'r--',lw=3,alpha=0.6)
    #plot( ones(2)* (mu - 2*sig), [-1e3,1e3],'r--',lw=3,alpha=0.6)
    #hist( cor[:,ang],bins=num_bin,normed=True,log=True)
    #show()


    total_del += len(bad1)
    print "total deleted shots:",total_del, 'remaining angles:', num_phi/2 - edge - ang

out = h5py.File(argv[2],'w')
out.create_dataset('cor',data = cor)
out.create_dataset('fnames',data = fnames3)
out.close()


