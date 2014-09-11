from OrbFitLib import Orbit
import math

def print_conic(c):
    print("{0:.6e} {1:.6e} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6e} {7:.6e}".
            format(c.rp,c.e,c.i,c.O,c.w,c.M0,c.t0,c.mu))

rp = 5.0
e  = 0.6
i  = 0.4
O  = 1.3
w  = 4.8
M0 = 3.4
t0 = 0
t1 = 100
mu = 100

e0 = [rp,e,i,O,w,M0,t0,mu]
c0 = Orbit()
c0.setup_elements(e0)
print_conic(c0)

r,v = c0.rv(t1)
c1 = Orbit()
c1.setup_rv(t1,mu,r,v)
print_conic(c1)

e2 = c1.elements(t0)
c0.setup_elements(e2)
print_conic(c0)

err = 0 
for diff in [e2[i]-e0[i] for i in range(len(e0))]: err += diff*diff
print(math.sqrt(err/8.))
