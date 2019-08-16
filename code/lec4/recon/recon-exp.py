from pylab import *
from numpy import *
import scipy.optimize

style.use('postgkyl.mplstyle')

# input values for f0 + x*f1

def func(an, f0, f1):
    a0 = an[0]; a1 = an[1]
    rhs0 = (exp(a1+a0) - exp(a0-a1))/a1
    rhs1 = ((a1-1)*exp(a0+a1) + (a1+1)*exp(a0-a1))/a1**2

    return rhs0-2*f0, rhs1-2.0/3.0*f1

def calcExpFit(f0, f1):
    # compute g0 and g1 for f0=1, f1=1.0 with initial guess 1.0, 0.01
    aout = scipy.optimize.fsolve(func, [1.0, 0.01], args=(f0, f1))
    g0 = aout[0]
    g1 = aout[1]

    return g0, g1

X = linspace(-1, 1, 100)

# fit exponential
f0 = 1.0
f1 = 1.0
g0, g1 = calcExpFit(f0, f1)
figure(1)
plot(X, f0+X*f1, 'r-')
plot(X, exp(g0+X*g1), 'k-')
gca().set_ylim([-1, 6])
grid()

savefig('exp-fit-1.png')

f0 = 1.0
f1 = 2.0
g0, g1 = calcExpFit(f0, f1)
figure(2)
plot(X, f0+X*f1, 'r-')
plot(X, exp(g0+X*g1), 'k-')
grid()
savefig('exp-fit-2.png')

show()

