from pylab import *
from healpy import *
from lognormal_mocks import *

seed(329968)


Nobj = 20000


nside = 512
npix = 12*nside**2

Nl = 768

Ntheta = 2048



rhobar = array([float(Nobj)/npix])

print(rhobar)

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
plot(ClA,label='input Cl')
plot(rhobar**2 * Clgauss[0,0],label='Gaussianized Cl')
plot(Clm,label='observed Cl')
loglog()

legend()


figure()
hist(m,200)

show()
