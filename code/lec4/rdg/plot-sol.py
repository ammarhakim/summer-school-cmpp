from pylab import *
import tables
import postgkyl as pg

style.use('postgkyl.mplstyle')

def exactSol(X, Y, t):
    return exp(-2*t)*sin(X)*cos(Y)

data = pg.GData("rdg-diffuse-2d_q_1.h5")
dg = pg.data.GInterpModal(data, 1, "ns")
XX, q0 = dg.interpolate()

Xhr = linspace(0, 2*pi, 101)
Yhr = linspace(0, 2*pi, 101)
XXhr, YYhr = meshgrid(Xhr, Yhr)
fhr = exactSol(XXhr, YYhr, 1.0)

figure(1)
pcolormesh(Xhr, Yhr, fhr)
colorbar()
figure(2)
pcolormesh(XX[0], XX[1], transpose(q0[:,:,0]))
colorbar()

show()
