## Solution of the 2D wave equation
from __future__ import division
import copy
from numpy import *
from pylab import *
from time import clock,time
from scipy import *
from mpl_toolkits.mplot3d import Axes3D

start = time()
def wave2d(c,dt,dx,dy,T,n):
	up = zeros((n,n))
	u = up.copy()
	um = up.copy()
	rx = c*(dt/dx)
	ry = c*(dt/dy)

	
	# Solution at t = 0
	start = time()
	for i in range(1,n):
		for j in range(1,n):
			u[i,j] = 2*sin((pi/4.0)*i*dx)*sin((pi/4.0)*j*dy)
	# solution at time t = -1
	for i in range(1,n-1):
		for j in range(1,n-1):
			um[i,j] = u[i,j]+ 0.5*rx*(u[i-1,j]-2*u[i,j]+u[i+1,j]) + 0.5*ry*(u[i,j-1]-2*u[i,j]+u[i,j+1])
	
	# BOundary values for um
			um[0,j] = 0
			um[i,0] = 0
	um[0,n-1] = 0
	um[n-1,0] =0		
	# Solution at subsequent time step
	
	t = 0
	while (t<=T):
		t+=0.5*dt
		for i in range(1,n-1):
			for j in range(1,n-1):
				up[i,j] = 2*u[i,j]-um[i,j] + 0.5*rx*(u[i-1,j]-2*u[i,j]+u[i+1,j]) + 0.5*ry*(u[i,j-1]-2*u[i,j]+u[i,j+1])
				
		#  Boundary condition and Update for next time level
				up[0,j] = 0
				up[i,0] = 0
		up[0,n-1] = 0
		up[n-1,0] =0
		um = u.copy(); u = up.copy()
	elapsed = time()-start
	print "Elapsed time is: ",elapsed,"seconds for n = ",n
	return up

c = 0.1
T=2

def exactSolution(x,y):
	ue = 2*sin((pi/4.0)*x)*sin((pi/4.0)*y)*cos((sqrt(2)*pi*c*T)/4)
	return ue

n = 1000
dx = 10.0/n
dy = 10.0/n
dt = 10./n

x = linspace(-5,5,n)
y = linspace(-5,5,n)
X,Y = meshgrid(x,y)

ue = exactSolution(X,Y)
ua = wave2d(c,dt,dx,dy,T,n)
#print "Elapsed time is:,",elapsed,"second"
f#ig=figure()
#ax = Axes3D(fig)
#ax.plot_surface(X,Y,ua)
#show()






				
			
			
		
