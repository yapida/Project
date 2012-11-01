## Solution of the 2D wave equation
from __future__ import division
import copy
from numpy import *
from pylab import *
from time import clock,time
from scipy import *
from mpl_toolkits.mplot3d import Axes3D

#start = time()
def wave2d(c,dt,dx,dy,T,n):
	up = zeros((n,n))
	u = up.copy()
	um = up.copy()
	rx = (c*(dt/dx))**2
	ry = (c*(dt/dy))**2

	
	# Solution at t = 0
	
	for i in range(1,n):
		for j in range(1,n):
			u[i,j] = 2*sin((pi/4.0)*i*dx)*sin((pi/4.0)*j*dy)
			#u[i][j] = sin(2*pi*(i*dx+j*dy))
	start = time()
	um[1:n-1,1:n-1] = u[1:n-1,1:n-1]+ 0.5*rx*(u[0:n-2,1:n-1]-2*u[1:n-1,1:n-1]+u[2:n,1:n-1]) + 0.5*ry*(u[1:n-1,0:n-2]-2*u[1:n-1,1:n-1]+u[1:n-1,2:n])
	# BOundary values for um
	um[0,1:n-2] = 0
	um[1:n-2,0] = 0
	um[0,n-1] = 0
	um[n-1,0] =0		
	# Solution at subsequent time step
	
	t = 0
	while (t<=T):
		t+=dt
		up[1:n-1,1:n-1]=2*u[1:n-1,1:n-1]-um[1:n-1,1:n-1]+0.5*rx*(u[0:n-2,1:n-1]-2*u[1:n-1,1:n-1]+u[2:n,1:n-1])+\
		0.5*ry*(u[1:n-1,0:n-2]-2*u[1:n-1,1:n-1]+u[1:n-1,2:n])
				
		#  Boundary condition and Update for next time level
		up[0,1:n-2] = 0
		up[1:n-2,0] = 0
		up[0,n-1] = 0
		up[n-1,0] =0
		um = u.copy(); u = up.copy()
	elapsed = time()-start
	print"Elapdes time", elapsed,"second for n = ",n
	return up
#elapsed = time()-start
c = 1.5
T=2

def exactSolution(x,y):
	ue = 2*sin((pi/4.0)*x)*sin((pi/4.0)*y)*cos((sqrt(2)*pi*c*T*0.25))
	#ve = sin(2*pi*(x+y))*cos(4*pi*T*c)
	return ue

n = 600
dx = 1./n
dy = 1./n
dt = 1./n

x = linspace(-1,1,n)
y = linspace(-1,1,n)
X,Y = meshgrid(x,y)

ue = exactSolution(X,Y)
ua = wave2d(c,dt,dx,dy,T,n)
#error = ue-ua
#print "Elapsed time :",elapsed,"second"
#fig=figure()
#ax = Axes3D(fig)
#ax.plot_surface(X,Y,error)
#show()






				
			
			
		
