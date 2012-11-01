/* ****************************************************************************************************
 * 2D WAVE-EQUATION using  explicite finite difference method
 *
 * 2D WAVE-EQUATION: u_{tt} = c *c* ( u_{xx} + u_{yy} )
 *
 * Testing for:
 *      Exact solution: u(t,x,y) =  2sin(0.25pix)(sin(0.25piy))cos(sqrt(2)pi*0.25*t*c)
 *      Boundary condition: u(t,0,y)=u(t,x,0)=0
 *
 *
 * *****************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifndef PI 
#define PI 4.*atan(1.) // 3.14159265358979323846
#endif
void Init0(double **u,double x,double y, double dx, double dt,int M,int N);
void Init1(double **u, double **um,double dx,double dt, double c, int M,int N);
void wave2d(double *u_buffer, double *um_buffer, double *up_buffer, double **u, double **um, double **up, double dt, double dx, double dy);
void deallocate(double *u_buffer,double *um_buffer,double *up_buffer,double **u,double **um,double **up);

int main(int argc, char *argv[]) 
{
  double c = 0.0001;
  int L = 1000; //Length of the interval
  int i, j; //Iterators 
  double t = 0.; //Current timestep
  double x = 0.; //Current grid point in x direction initialized
  double y = 0.; //Current grid point in y direction
  double dt; //The size of the time-steps
  double dx; //The size of the the spatial intervals (spacing between the x_i and x_(i+1))
  double dy; //The size of the the spatial intervals (spacing between the y_j and y_(j+1))
  double *u_buffer, *um_buffer, *up_buffer; /* The 1D continous allocated memory for the 3 2D-grids of u(t) ,um(t-1) and up(t+1) */
  double **u, **um, **up; //Pointers to the different rows in the 3 2D grids for a given timestep
  int M = L-1;
  int N = L-1;
  int T = 1000;
  dx = L/(M+1.);
  dy = L/(N+1.);
  dt = T/L;
  clock_t start,end;
  double elapsed;
  
 
  /* Allocating of 3D arrays up, u, um */
  if (( um_buffer = (double *) malloc( sizeof(double) * (M+2)*(N+2) ) ) == NULL ) {
    fprintf(stderr, "Out of memory\n");
    exit(0);
  }
  if (( u_buffer = (double *) malloc( sizeof(double) * (M+2)*(N+2) ) ) == NULL ) {
    fprintf(stderr, "Out of memory\n");
    exit(0);
  }
  if (( up_buffer = (double *) malloc( sizeof(double) * (M+2)*(N+2) ) ) == NULL ) {
    fprintf(stderr, "Out of memory\n");
    exit(0);
 }  
  /**/
  if (( um = (double **) malloc( sizeof(double) * (N+2) ) ) == NULL ) {
    fprintf(stderr, "Out of memory\n");
    exit(0);
  }
  for (i=0; i<N+2; i++) {
    um[i] = &(um_buffer[i*(M+2)]);
  }
  if (( u = (double **) malloc( sizeof(double) * (N+2) ) ) == NULL ) {
    fprintf(stderr, "Out of memory\n");
    exit(0);
  }
  for (i=0; i<N+2; i++) {
    u[i] = &(u_buffer[i*(M+2)]);
  }
  if (( up = (double **) malloc( sizeof(double) * (N+2) ) ) == NULL ) {
    fprintf(stderr, "Out of memory\n");
    exit(0);
  }
  for (i=0; i<N+2; i++) {
    up[i] = &(up_buffer[i*(M+2)]);
  }


   start = clock();
   // Initialize for outer points
   Init0(u,x,y,dx,dt,M,N);
   
   // Initialize for inner point
   Init1(u,um,dx,dt,c,M,N);
   
   // Solve the wave eqation
   wave2d(u_buffer,um_buffer,up_buffer, u,um,up, dt, dx,dy);
   end = clock();
   elapsed = (double)(end-start)/CLOCKS_PER_SEC;
   
   printf("Elapsed time is: %f\n seconds ", elapsed );
   // Deallocate memory
   deallocate(u_buffer,um_buffer,up_buffer,u,um,up);
   
   return 0;
}
  



   
  /* Enforcing initial condition 1*/
void Init0(double **u,double x,double y,double dx,double dy,int M, int N )
{  
  int i,j;
  for (j=0; j<=N+1; j++)
  {
	y = j*dy;
	for (i=0; i<=M+1; i++)
	{
		x = i*dx;
		//u[j][i] = sin(2*PI*(x+y));
		u[j][i] = 2*sin(0.25*PI*x)*sin(0.25*PI*y);
	}
  }
}




  /* Initial condition for the interior points */
void Init1(double **u, double **um,double dx,double dt, double c, int M, int N)
{  
  int i,j;
  double r =(c*(dt/dx))*((c*(dt/dx)));
  
  for (j=1; j<=N; j++)
  {
	um[j][0]=0;
	for (i=1; i<=M; i++)
	{
		um[0][i]=0;
		um[j][i] = u[j][i]\
			 +0.5*r*(u[j][i-1]-2*u[j][i]+u[j][i+1])\
			 +0.5*r*(u[j-1][i]-2*u[j][i]+u[j+1][i]);
	}
	
  }

}


void wave2d(double *u_buffer, double *um_buffer, double *up_buffer, double **u, double **um, double **up, double dt, double dx, double dy)
{
   int i,j;
   double **tmp;// Temporary pointer
   double v;
   double L2_norm;
   double inner_P;
   
   double c = 0.0001;
   
   double r1 = c*c*(dt*dt)/(dx*dx);
   double y=0.;
   double x=0.;
   int L =1000;
   int M = L-1;
   int N = L-1;
   int T = 1000;
   dx = L/(M+1.);
   dy = L/(N+1.);
  // double t_max = 10.;
   dt = T/L;
   

  /* Main time loop */
  int 
   t=0;
   while(t<T)
  { // Start while
	t += dt;
	for (i=0; i<=M+1; i++)
	{
		x=i*dx;
		up[0][i]=0; // BC at y=0
	}

	for(j=1; j<=N; j++)
	{
		y=j*dy;
		up[j][0]=0; // BC at x=0
		for(i=1; i<=M; i++)
		{
			up[j][i] = 2*u[j][i]-um[j][i]\
			         +0.5*r1*(u[j-1][i]-2*u[j][i]+u[j+1][i])\
			         +0.5*r1*(u[j][i-1]-2*u[j][i]+u[j][i+1]);	
		}

	}
	
	/* Update data for next time step */
	tmp = um;
	um = u;
	u = up;
	up = tmp;
  }// End while
	
	// Print solution
       for (j=0; j<=M+1; j++)
       {
   	for (i=0; i<=N+1; i++)
   		printf("u(%d,%d) = %g\n" ,i,j,up[j][i]);
       }
       
	// Print L2 norm
	c=0.0001;
	 for(j=1; j<=N; j++)
       {
		y=j*dy;
		for(i=1;i<=M; i++)
		{
			x=i*dx;
			 //v=sin(2*PI*(x+y))*cos(4*PI*T*c);
			 v = 2*sin(0.25*PI*x)*sin(0.25*PI*y)*cos(sqrt(2)*c*T*0.25*PI);
			 inner_P += (v-up[j][i])*(v-up[j][i]);
		}
		
       }
       L2_norm=sqrt(inner_P/(M+N));
       printf("L2_norm:%g\n",L2_norm);



   
}// End function wave2d


// Memory deallocation function
void deallocate(double *u_buffer,double *um_buffer,double *up_buffer,double **u,double **um,double **up)
{
   free(u_buffer);
   free(um_buffer);
   free(up_buffer);
   free(u);
   free(um);
   free(up);	
}


