import sys

from pylab import *

from healpy import *

sys.path.reverse()
sys.path.append('build/lib.linux-x86_64-2.7/')
sys.path.reverse()
from lognormal_mocks import *

helloworld()

seed(329968)


Nobj = 20000


nside = 512
npix = 12*nside**2

Nl = 768

Ntheta = 2048



rhobar = array([float(Nobj)/npix])

print rhobar

l = arange(0.,Nl)

#ClA = (exp(-l**2/2./20.**2) + 0.3 * exp(-(l-50)**2/2./20.**2))*(l+1)**-2
ClA = (l**0.5 * exp(-l**2/2./320.**2))*1e-10


Cl = array([[ ClA ]])


gaussbar, Clgauss = lognormal_mocks_stats(rhobar,Cl,Ntheta)


mgauss = synfast(Clgauss[0,0],nside)

m = exp(gaussbar[0] + mgauss)

Clm = anafast(m)

close('all')

mollview(m)

figure()
plot(ClA)
plot(rhobar**2 * Clgauss[0,0])
plot(Clm)
#ylim(ymax=amax(Cl)*1.1)
loglog()
#xlim(0,Nl)


figure()
hist(m,200)

show()
