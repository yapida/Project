## Solving the 1D wave equation u_tt = c*c(u_xx + u_yy),
## u(x,0) = sin(2pix)
## u(0,t) = u(1,t) = 0, u_t(x,0) = g =0,
import copy
from numpy import *
from pylab import *
import pylab
from time import clock,time
from scipy import *
from mpl_toolkits.mplot3d import Axes3D
def wave1d(c,dt,dx,Tstop,n):
	r = c*(dt/dx)
	
## setting boundary conditions and initial conditions
## initial condition u(x,0) = sin(2pix)
## solution at t = 0 for all space points
	up = zeros(n)  # Solution u(t+1)
	u = up.copy()  # Solution u(t)
	um = up.copy() # Solution u(t-1)
	# Solution at t=0
	start = time()
	for j in range(1,n):
		u[j] = sin(2*pi*j*dx)
	for j in range(1,n-1):
		um[0] = 0
		um[n-1] = 0
		um[j] = u[j]  + 0.5*r**2*(u[j-1] - 2*u[j] + u[j+1])
	
	
## solution at subsequent time step, t 
	t = 0
	while (t <= T):
		t+=dt
		for j in range(1,n-1):
			up[0] = 0
			up[n-1] = 0
			up[j] = 2*u[j] - um[j] + r**2*(u[j+1]-2*u[j] + u[j-1])
		# Update for next time level
		um = u.copy(); u = up.copy()
	
	elapsed = time()-start
	print "Elapsed time:",elapsed,"second, for n = ",n
				
	return up


#Exact solution
T = 2
c = 1.5
def exactSolution(x):
	return 0.5*( sin(2*pi*(x-c*T)) + sin(2*pi*(x+c*T)) )

n=900
dt = 1.0/n
dx = 1.0/n

x = linspace(0,1,n)
#T = linspace(0,1,n)
#X,Y = meshgrid(x,T)

ue = exactSolution(x)
ua = wave1d(c,dt,dx,T,n)
#ua.resize([n,n])
#fig = figure()
#ax = Axes3D(fig)

#ax.plot_surface(X,Y,ua)
#pylab.plot(x,ua,x,ue)

#pylab.show()
#print ue, clock()



			
			

		
	

	
	

	

