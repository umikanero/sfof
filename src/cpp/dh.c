/* ----------------------------------------------------------------------------

   Cosmography routines in C.

-------------------------------------------------------------------------------

   Use at your own risk.  Note that, where necessary, integration is
   done with stepsize roughly MINSTEP and that the integrator is not
   pretty.

-------------------------------------------------------------------------------
   David W. Hogg
   hogg@ias.edu
   Modificado por Eduardo Cypriano (03/02/2000)
---------------------------------------------------------------------------- */
#define NR_END 0 
#define FREE_ARG char*


double angdidis(double,double,double);
double lumdis(double,double,double);
void add_protocluster_to_cluster (struct protocluster *, struct cluster *);
void clean_repetitions (struct cluster *);
void convert_decimal_to_hms (double , int *, int *, float *);
void convert_decimal_to_dms (double , char *, int *, int *, float *);

double angdidis(double z,double OmegaM,double OmegaL);
double angdidis2(double z1,double z2,double OmegaM,double OmegaL);
double comdis(double z,double OmegaM,double OmegaL);
double comvol(double z,double OmegaM,double OmegaL);
double dcomdisdz(double z,double OmegaM,double OmegaL);
double dlookbackdz(double z,double OmegaM,double OmegaL);
double doptdepthdz(double z,double OmegaM,double OmegaL);
double dpropmotdisdz(double z,double OmegaM,double OmegaL);
double intcomvol(double z,double OmegaM,double OmegaL);
double lookback(double z,double OmegaM,double OmegaL);
double optdepth(double z,double OmegaM,double OmegaL);
double propmotdis(double z,double OmegaM,double OmegaL);

#include <math.h>
#include<stdio.h>
#define PI       (3.141592653589793238462643)
#define TINY     (1.0e-16)
#define MINSTEP  (0.01)


/* ----------------------------------------------------------------------------
   angdidis
-------------------------------------------------------------------------------
   This function calculates the angular diameter distance d_A as a
   function of z, Omega_M and Omega_L in a matter-dominated universe,
   using the function propmotdis().  H0=c=1.
---------------------------------------------------------------------------- */
double angdidis(double z,double OmegaM,double OmegaL)
{
  return propmotdis(z,OmegaM,OmegaL)/(1.0+z) ;
}


/* ----------------------------------------------------------------------------
   angdidis2
-------------------------------------------------------------------------------
   This function calculates the angular diameter distance d_A from z1
   to z2 as a function of Omega_M and Omega_L in a matter-dominated
   universe, using the function propmotdis().  H0=c=1.
---------------------------------------------------------------------------- */
double angdidis2(double z1,double z2,double OmegaM,double OmegaL)
{
  double y1,y2,y12,OmegaR ;

  OmegaR= 1.0-OmegaM-OmegaL ;
  if(OmegaR<-TINY) printf("WARNING: angdidis2() does not work at Omega_R<0!") ;
  y1= propmotdis(z1,OmegaM,OmegaL) ;
  y2= propmotdis(z2,OmegaM,OmegaL) ;
  y12= y2*sqrt(1.0+y1*y1*OmegaR)-y1*sqrt(1.0+y2*y2*OmegaR) ;

  return y12/(1.0+z2) ;
}


/* ----------------------------------------------------------------------------
   comdis
-------------------------------------------------------------------------------
   This function calculates the line-of-sight comoving distance d_C as
   a function of z, Omega_M and Omega_L in a matter-dominated
   universe, using dcomdisdz().  H0=c=1.
---------------------------------------------------------------------------- */
double comdis(double z,double OmegaM,double OmegaL)
{
  int nsteps ;
  double dz,dC,zz ;

  nsteps= ((int) (z/MINSTEP))+1 ;
  dz= z/((double) nsteps) ;
  dC= 0.0 ;
  for(zz=0.5*dz; zz<z; zz+=dz) dC += dz*dcomdisdz(zz,OmegaM,OmegaL) ;

  return dC ;
}


/* ----------------------------------------------------------------------------
   comvol
-------------------------------------------------------------------------------
   This function calculates the all-sky comoving volume V as a
   function of z, Omega_M and Omega_L in a matter-dominated universe.
   Formulae from Carrol, Press & Turner, 1992, and my own calculation.
---------------------------------------------------------------------------- */
double comvol(double z,double OmegaM,double OmegaL)
{
  double V,dM,OmegaK,sqrtOmegaK ;

  OmegaK= 1.0-OmegaM-OmegaL ;
  sqrtOmegaK= sqrt(fabs(OmegaK)) ;
  dM= propmotdis(z,OmegaM,OmegaL) ;

  if(OmegaK < -TINY)
    V= (dM*sqrt(1.0+OmegaK*dM*dM)-asin(dM*sqrtOmegaK)/sqrtOmegaK)
      /(2.0*OmegaK) ;
  else if(OmegaK > TINY)
    V= (dM*sqrt(1.0+OmegaK*dM*dM)-asinh(dM*sqrtOmegaK)/sqrtOmegaK)
      /(2.0*OmegaK) ;
  else
    V= dM*dM*dM/3.0 ;

  return 4.0*PI*V ;
}


/* ----------------------------------------------------------------------------
   dcomdisdz
-------------------------------------------------------------------------------
   This function calculates the differential line-of-sight comoving
   distance dD_c/dz as a function of z, Omega_M and Omega_L in a
   matter-dominated universe.  H0=c=1.
---------------------------------------------------------------------------- */
double dcomdisdz(double z,double OmegaM,double OmegaL)
{
  return (1.0/sqrt((1.0+z)*(1.0+z)*(1.0+OmegaM*z)-z*(2.0+z)*OmegaL)) ;
}


/* ----------------------------------------------------------------------------
   dcomvoldz
-------------------------------------------------------------------------------
   This function calculates the one-steradian differential comoving
   volume dV/dz as a function of z, Omega_M and Omega_L in a
   matter-dominated universe.  Formulae from Carrol, Press & Turner,
   1992, Kolb & Turner, 1990, and my own calculation.  H0=c=1.
---------------------------------------------------------------------------- */
double dcomvoldz(double z,double OmegaM,double OmegaL)
{
  double dM,OmegaK,ddMdz ;

  OmegaK= 1.0-OmegaM-OmegaL ;
  dM= propmotdis(z,OmegaM,OmegaL) ;
  ddMdz= dpropmotdisdz(z,OmegaM,OmegaL) ;

  return dM*dM*ddMdz/sqrt(1.0+OmegaK*dM*dM) ;
}


/* ----------------------------------------------------------------------------
   dlookbackdz
-------------------------------------------------------------------------------
   This function calculates the change in lookback time dt/dz with
   redshift z as a function of z, Omega_M and Omega_L in a
   matter-dominated universe.  Formula from Carrol, Press & Turner,
   1992.  H0=c=1.
---------------------------------------------------------------------------- */
double dlookbackdz(double z,double OmegaM,double OmegaL)
{
  return 1.0/((1.0+z)*sqrt((1.0+z)*(1.0+z)*(1.0+OmegaM*z)-z*(2.0+z)*OmegaL)) ;
}


/* ----------------------------------------------------------------------------
   doptdepthdz
-------------------------------------------------------------------------------
   This function calculates the change in optical depth dtau/dz with
   redshift z as a function of z, Omega_M and Omega_L in a
   matter-dominated universe.  Formula from Peebles, 1993.
   H0=c=sigma=n=1.
---------------------------------------------------------------------------- */
double doptdepthdz(double z,double OmegaM,double OmegaL)
{
  return (1.0+z)*(1.0+z)/
    sqrt((1.0+z)*(1.0+z)*(1.0+OmegaM*z)-z*(2.0+z)*OmegaL) ;
}


/* ----------------------------------------------------------------------------
   dpropmotdisdz
-------------------------------------------------------------------------------
   This function calculates the derivative of the proper motion
   distance d_M with respect to redshift z as a function of z, Omega_M
   and Omega_L in a matter-dominated universe.  Formula from Carrol,
   Press & Turner, 1992.  This function also requires the function
   propmotdis(). H0=c=1.
---------------------------------------------------------------------------- */
double dpropmotdisdz(double z,double OmegaM,double OmegaL)
{
  double ddMdz,OmegaK,dM;

  ddMdz = 1.0/sqrt((1.0+z)*(1.0+z)*(1.0+OmegaM*z)-z*(2.0+z)*OmegaL) ;

  OmegaK= 1.0-OmegaM-OmegaL ;
  if(OmegaK < -TINY){
    dM= propmotdis(z,OmegaM,OmegaL) ;
    ddMdz= sqrt(1.0-OmegaK*dM*dM)*ddMdz ;
  }else if(OmegaK > TINY){
    dM= propmotdis(z,OmegaM,OmegaL) ;
    ddMdz= sqrt(1.0+OmegaK*dM*dM)*ddMdz ;
  }

  return ddMdz ;
}


/* ----------------------------------------------------------------------------
   intcomvol
-------------------------------------------------------------------------------
   This function calculates the all-sky comoving volume V as a
   function of z, Omega_M and Omega_L in a matter-dominated universe
   by integrating dcomvoldz().  It was written to test comvol().
---------------------------------------------------------------------------- */
double intcomvol(double z,double OmegaM,double OmegaL)
{
  int nsteps ;
  double dz,zz,V;

  nsteps= ((int) (z/MINSTEP))+1 ;
  dz= z/((double) nsteps) ;
  V= 0.0 ;
  for(zz=0.5*dz;zz<z;zz+=dz) V+= dz*dcomvoldz(zz,OmegaM,OmegaL) ;

  return 4.0*PI*V ;
}


/* ----------------------------------------------------------------------------
   lookback
-------------------------------------------------------------------------------
   This function calculates the lookback time t(0)-t(z) as a function
   of z, OmegaM and OmegaL by integrating the output of dlookbackdz.
---------------------------------------------------------------------------- */
double lookback(double z,double OmegaM,double OmegaL)
{
  int nsteps ;
  double t,zz,dz ;
  nsteps= ((int) (z/MINSTEP))+1 ;
  dz= z/((double) nsteps) ;
  t= 0.0 ;
  for(zz=0.5*dz; zz<z; zz+=dz) t+= dlookbackdz(zz,OmegaM,OmegaL)*dz ;
  return t ;
}


/* ----------------------------------------------------------------------------
   lumdis
-------------------------------------------------------------------------------
   This function calculates the luminosity distance d_L as a function
   of z, Omega_M and Omega_L in a matter-dominated universe, using the
   function propmotdis().  H0=c=1.
---------------------------------------------------------------------------- */
double lumdis(double z,double OmegaM,double OmegaL)
{
  return propmotdis(z,OmegaM,OmegaL)*(1.0+z) ;
}


/* ----------------------------------------------------------------------------
   optdepth
-------------------------------------------------------------------------------
   This function calculates the optical depth tau as a function
   of z, OmegaM and OmegaL by integrating the output of doptdepthdz.
   Again, H0=c=sigma=n=1.
---------------------------------------------------------------------------- */
double optdepth(double z,double OmegaM,double OmegaL)
{
  int nsteps ;
  double tau,zz,dz ;

  nsteps= ((int) (z/MINSTEP))+1 ;
  dz= z/((double) nsteps) ;
  tau= 0.0 ;
  for(zz=0.5*dz; zz<z; zz+=dz) tau+= doptdepthdz(zz,OmegaM,OmegaL)*dz ;
  return tau ;
}


/* ----------------------------------------------------------------------------
   propmotdis
-------------------------------------------------------------------------------
   This function calculates the proper motion distance d_M as a
   function of z, Omega_M and Omega_L in a matter-dominated universe.
   Formulae from Carrol, Press & Turner, 1992, Kolb \& Turner, 1990,
   and my own derivation.  Makes use of comdis().  H0=c=1.
---------------------------------------------------------------------------- */
double propmotdis(double z,double OmegaM,double OmegaL)
{
  double dM,q0,OmegaK,sqrtOmegaK;

  if(fabs(OmegaM)<TINY && fabs(OmegaL)<TINY){
    dM= (z+0.5*z*z)/(1.0+z) ;

  }else if(fabs(OmegaL)<TINY){
    q0= 0.5*OmegaM ;
    dM= (z*q0+(q0-1.0)*(sqrt(2.0*q0*z+1.0)-1.0))/(q0*q0*(1.0+z)) ;

  }else{
    dM= comdis(z,OmegaM,OmegaL) ;
    OmegaK= 1.0-OmegaM-OmegaL ;
    sqrtOmegaK= sqrt(fabs(OmegaK)) ;
    if(OmegaK < -TINY) dM= sin(sqrtOmegaK*dM)/sqrtOmegaK ;
    else if(OmegaK > TINY) dM= sinh(sqrtOmegaK*dM)/sqrtOmegaK ;
  }

  return dM ;
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *fvector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror(const_cast<char *>("allocation failure in vector()"));
	return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


float splint(float xa[], float ya[], float y2a[], int n, float x)
{
        int klo,khi,k;
        float h,b,a,y;
	
	klo=1;
        khi=n;
        while (khi-klo > 1) {
                k=(khi+klo) >> 1;
                if (xa[k] > x) khi=k;
                else klo=k;
        }
        h=xa[khi]-xa[klo];
        if (h == 0.0) nrerror(const_cast<char *>("Bad xa input to routine splint"));
        a=(xa[khi]-x)/h;
        b=(x-xa[klo])/h;
        y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
	return y;
}

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[])
{
        int i,k;
        float p,qn,sig,un,*u;
                                                                                                               
        u=fvector(1,n-1);
        if (yp1 > 0.99e30)
                y2[1]=u[1]=0.0;
        else {
                y2[1] = -0.5;
                u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
        }
        for (i=2;i<=n-1;i++) {
                sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
                p=sig*y2[i-1]+2.0;
                y2[i]=(sig-1.0)/p;
                u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
                u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
        }
        if (ypn > 0.99e30)
                qn=un=0.0;
        else {
                qn=0.5;
                un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
        }
        y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
        for (k=n-1;k>=1;k--)
                y2[k]=y2[k]*y2[k+1]+u[k];
        free_vector(u,1,n-1);
}
