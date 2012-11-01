## Solving the 1D wave equation u_tt = c*c(u_xx + u_yy)
## u(x,0)=sin(2pix)
## u(0,t) = u(1,t) = 0, u_t(x,0) = g =0,
import copy
from numpy import *
from pylab import *
import pylab
from time import clock, time
from scipy import *
from mpl_toolkits.mplot3d import Axes3D
def wave1d(c,dt,dx,T,n):
	r = c*(dt/dx)
	
## setting boundary conditions and initial conditions
## initial condition u(x,0) = sin(2pix)
## solution at t = 0 for all space points
	up = zeros(n)
	u = up.copy()
	um = up.copy()
	# Solution at t=0
	start = time()
	for j in range(1,n):
		u[j] = sin(2*pi*j*dx)
	# Solution a t = -1, Vecorized version
	um[0]=0
	um[n-1]=0
	um[1:n-1] = u[1:n-1]+0.5*r**2*(u[0:n-2]-2*u[1:n-1]+u[2:n])
## solution at subsequent time step; vectorized version
	t = 0
	while t<=T:
		t+=dt
		up[0] = 0
	 	up[n-1] = 0
		up[1:n-1] = 2*u[1:n-1] - um[1:n-1] + r**2*(u[0:n-2]-2*u[1:n-1] +u[2:n])
		# Update for next time level
		um = u.copy(); u = up.copy()
	
	
	elapsed = time()-start
	print "Elapsed time is:",elapsed,"seconds for n = ",n	
	return up


#Exact solution
c = 1.01
T = 1
def exactSolution(x):
	return 0.5*( sin(2*pi*(x-c*T)) + sin(2*pi*(x+c*T)) )



n=1000

dt = 5.0/n
dx = 5./n
x = linspace(0,5,n)

ue = exactSolution(x)
ua = wave1d(c,dt,dx,T,n)
pylab.plot(x,ua,x,ue)
pylab.show()

			
			

		
	

	
	

	

