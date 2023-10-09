
import pylab
import math
from rk import RK4

xdot = lambda t, y, z: math.cos(t)
print xdot
ydot = lambda t, y, z: -y
print ydot

rk4 = RK4(xdot, ydot)
t, y = rk4.solve([0, 1], .01, 2*math.pi)

pylab.plot(y[0], y[1])
pylab.show()
