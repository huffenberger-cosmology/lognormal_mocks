from pylab import *
from healpy import *
from lognormal_mocks import *
import numpy as np

close('all')

seed(329968)
nside = 1024
npix = 12*nside**2
omegapix = 4*pi/npix
Nl = 600

l = arange(40,600)

#start with some cl's
cltt = 1 * (l/80)**(-1+2)
clte = 0*l #2.77 * 315.4 * (l/80)**(-2.5+2)
clee = 315.4 * (l/80)**(-2.54+2)
cltb = 0*l # 1 * (l/80)**(-1+2)
cleb = 0*l # 1 * (l/80)**(-1+2)
clbb = 0.53 * 315.4 * (l/80)**(-2.54+2)

#clte = 2.77* 315.4 * (l/80)**(-2.5+2) (l from 40 to 600? LR71)
#clee = 315.4 * (l/80)**(-2.42+2)
#clbb = 0.53 * 315.4 * (l/80)**(-2.54+2)
#beta tt of dust is 1.48 (mean?)
# dust spectral index: btt = 1.5, bee = bbb = 1.59

# Truncate each spectrum to the same length

#these would be means, we can work with the wrong values first, and then adjust. (One of the means should be zero? but we work with a number first
# and then subtract
tbar = 100e6
ebar = 4000
bbar = 300e6

#kappabar = 100.0 # a large, fake mean kappa value makes the kappa field more gaussian, it is subtracted away in the end
rhobar = array([tbar,ebar,bbar])  # array of mean must be size Nmap

# Note the ngal factors make the statistics for the galaxy density field
Cl = array([[ cltt, clte, cltb ],
            [ clte, clee, cleb ],
            [ cltb, cleb, clbb ]])  # array of input Cl must be size (Nmap, Nmap, Nl).  If this is the wrong shape it will fail

Ntheta = 10000 # accuracy parameter

gaussbar, Clgauss = lognormal_mocks_stats(rhobar,Cl,Ntheta) # get the stats for the Gaussianized fields

almgauss = synalm( (Clgauss[0,0],Clgauss[1,1],Clgauss[2,2], Clgauss[0,1], Clgauss[1,2], Clgauss[0,2]),
               new=True) # make the correlated, gaussianized alms

t_gauss = alm2map(almgauss[0],nside) # synthesize the gaussianize maps
e_gauss = alm2map(almgauss[1],nside)
b_gauss = alm2map(almgauss[2], nside)

# exponentiate to get the non-gaussian maps
t = exp(gaussbar[0] + t_gauss) # This is the galaxy density field
e = exp(gaussbar[1] + e_gauss) - ebar
b = exp(gaussbar[2] + b_gauss) - bbar # This is the mean-zero kappa map
#subtract the mean for B

#STOPPED MODIFYING HERE

# evaluate the power spectra of the realized maps
Clt = anafast(t)
Cle = anafast(e)
Clb = anafast(b)

#might need anafast or synfast or alm thing for three cls or three maps

# Visualize the results
