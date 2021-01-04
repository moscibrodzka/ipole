#include "decs.h"

/*

  first geometry used in the model and then 
	HARM model specification routines 

	model-dependent geometry routines:
	        risco_calc
		rhorizon_calc
		gcov
		get_connection
*/

//Kerr metric
double risco_calc( int do_prograde )
{
  double Z1,Z2,sign,term1,term2 ;
  sign = (do_prograde) ? 1. : -1. ;
  term1 = pow(1. + a,1./3.);
  term2 = pow(1. - a,1./3.);
  Z1 = 1. + term1*term2*(term1 + term2);
  Z2 = sqrt(3.*a*a + Z1*Z1) ;
  return( 3. + Z2-sign*sqrt((3. - Z1)*(3. + Z1 + 2.*Z2))  );

}

//Kerr metric
double rhorizon_calc(int pos_sign)
{
  double sign;
  sign = (pos_sign) ? 1. : -1.;
  return(1. + sign*sqrt((1.-a)*(1.+a)) );
}

//function from IL group, coefficents needed to set up metric corrections for modified coodtinates 
#define MUNULOOP for(int mu=0;mu<NDIM;mu++) for(int nu=0;nu<NDIM;nu++)
void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{
    // Jacobian with respect to KS basis where X is given in
    // non-KS basis
    MUNULOOP dxdX[mu][nu] = 0.;

    //traditional coordinates with cut and with hslope
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][2] = M_PI + hslope*2*M_PI*cos(2*M_PI*X[2]); 
    dxdX[3][3] = 1.;

}


/* 
   Current metric: modified Kerr-Schild, squashed in theta
   to give higher resolution at the equator 
*/

//second part of the function is from IL group
void gcov_func(double *X, double gcov[][NDIM])
{

  double sth, cth, s2, rho2;
  double r, th;
  double gcov_ks[NDIM][NDIM];
  MUNULOOP gcov_ks[mu][nu] = 0.;

  bl_coord(X,&r,&th);

  cth = cos(th);
  sth = sin(th);
  s2 = sth * sth;
  rho2 = r * r + a * a * cth * cth;

  //Kerr-Schield metric tensor
  gcov_ks[0][0] = -1. + 2. * r / rho2 ;
  gcov_ks[0][1] = 2. * r / rho2 ;
  gcov_ks[0][3] = -2. * a * r * s2 / rho2;
  gcov_ks[1][0] = gcov_ks[0][1];
  gcov_ks[1][1] = 1. + 2. * r / rho2 ;
  gcov_ks[1][3] = -a * s2 * (1. + 2. * r / rho2);
  gcov_ks[2][2] = rho2 ;
  gcov_ks[3][0] = gcov_ks[0][3];
  gcov_ks[3][1] = gcov_ks[1][3];
  gcov_ks[3][3] = s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2));
    
  // convert from ks metric to a modified one using Jacobian
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);
  MUNULOOP
    {
      gcov[mu][nu] = 0;
      for (int lam = 0; lam < NDIM; ++lam) {
	for (int kap = 0; kap < NDIM; ++kap) {
	  gcov[mu][nu] += gcov_ks[lam][kap] * dxdX[lam][mu] * dxdX[kap][nu];
	}
      }
      
      if( isnan(gcov[mu][nu])){
	fprintf(stdout," gcov is nan here: r=%g th=%g X1=%g X2=%g -->>  %g %g \n",r,th,X[1],X[2],dxdX[2][1],dxdX[2][2]);
	exit(1);
      }
      
    }
  
}
#undef MUNULOOP




/* 

   connection calculated analytically for modified Kerr-Schild
	coordinates 
   this gives the connection coefficient
	\Gamma^{i}_{j,k} = conn[..][i][j][k]
   where i = {1,2,3,4} corresponds to, e.g., {t,ln(r),theta,phi}

*/


void get_connection(double X[NDIM], double lconn[NDIM][NDIM][NDIM])
{
  
    double r1,r2,r3,r4,sx,cx;
    double th,dthdx2,dthdx22,d2thdx22,sth,cth,sth2,cth2,sth4,cth4,s2th,c2th;
    double a2,a3,a4,rho2,irho2,rho22,irho22,rho23,irho23,irho23_dthdx2;
    double fac1,fac1_rho23,fac2,fac3,a2cth2,a2sth2,r1sth2,a4cth4;
	
    r1 = exp(X[1]);
    r2 = r1*r1;
    r3 = r2*r1;
    r4 = r3*r1;
  
    sx = sin(2.*M_PI*X[2]);
    cx = cos(2.*M_PI*X[2]);
  
    /* HARM-3D MKS */
    th =  M_PI*X[2] + hslope*sx;
    dthdx2 = M_PI + 2.*M_PI*hslope*cx;
    d2thdx22 = -4.*M_PI*M_PI*hslope*sx;	/* d^2(th)/dx2^2 */
    
  
    dthdx22 = dthdx2*dthdx2;
	
    sth = sin(th);
    cth = cos(th);
    sth2 = sth*sth;
    r1sth2 = r1*sth2;
    sth4 = sth2*sth2;
    cth2 = cth*cth;
    cth4 = cth2*cth2;
    s2th = 2.*sth*cth;
    c2th = 2*cth2 - 1.;
    
    a2 = a*a;
    a2sth2 = a2*sth2;
    a2cth2 = a2*cth2;
    a3 = a2*a;
    a4 = a3*a;
    a4cth4 = a4*cth4;
  
    rho2 = r2 + a2cth2;
    rho22 = rho2*rho2;
    rho23 = rho22*rho2;
    irho2 = 1./rho2;
    irho22 = irho2*irho2;
    irho23 = irho22*irho2;
    irho23_dthdx2 = irho23/dthdx2;
    
    fac1 = r2 - a2cth2;
    fac1_rho23 = fac1*irho23;
    fac2 = a2 + 2*r2 + a2*c2th;
    fac3 = a2 + r1*(-2. + r1);
    
    lconn[0][0][0] = 2.*r1*fac1_rho23;
    lconn[0][0][1] = r1*(2.*r1+rho2)*fac1_rho23;
    lconn[0][0][2] = -a2*r1*s2th*dthdx2*irho22;
    lconn[0][0][3] = -2.*a*r1sth2*fac1_rho23;
    
    lconn[0][1][0] = lconn[0][0][1];
    lconn[0][1][1] = 2.*r2*(r4 + r1*fac1 - a4cth4)*irho23;
    lconn[0][1][2] = -a2*r2*s2th*dthdx2*irho22;
    lconn[0][1][3] = a*r1*(-r1*(r3 + 2*fac1) + a4cth4)*sth2*irho23;
  
    lconn[0][2][0] = lconn[0][0][2];
    lconn[0][2][1] = lconn[0][1][2];
    lconn[0][2][2] = -2.*r2*dthdx22*irho2;
    lconn[0][2][3] = a3*r1sth2*s2th*dthdx2*irho22;

    lconn[0][3][0] = lconn[0][0][3];
    lconn[0][3][1] = lconn[0][1][3];
    lconn[0][3][2] = lconn[0][2][3];
    lconn[0][3][3] = 2.*r1sth2*(-r1*rho22 + a2sth2*fac1)*irho23;
    
    lconn[1][0][0] = fac3*fac1/(r1*rho23);
    lconn[1][0][1] = fac1*(-2.*r1 + a2sth2)*irho23;
    lconn[1][0][2] = 0.;
    lconn[1][0][3] = -a*sth2*fac3*fac1/(r1*rho23);

    lconn[1][1][0] = lconn[1][0][1];
    lconn[1][1][1] = (r4*(-2. + r1)*(1. + r1) + a2*(a2*r1*(1. + 3.*r1)*cth4 + a4cth4*cth2 + r3*sth2 + r1*cth2*(2.*r1 + 3.*r3 - a2sth2)))*irho23;
    lconn[1][1][2] = -a2*dthdx2*s2th/fac2;
    lconn[1][1][3] = a*sth2*(a4*r1*cth4 + r2*(2*r1 + r3 - a2sth2) + a2cth2*(2.*r1*(-1. + r2) + a2sth2))*irho23;
    
    lconn[1][2][0] = lconn[1][0][2];
    lconn[1][2][1] = lconn[1][1][2];
    lconn[1][2][2] = -fac3*dthdx22*irho2;
    lconn[1][2][3] = 0.;

    lconn[1][3][0] = lconn[1][0][3];
    lconn[1][3][1] = lconn[1][1][3];
    lconn[1][3][2] = lconn[1][2][3];
    lconn[1][3][3] = -fac3*sth2*(r1*rho22 - a2*fac1*sth2)/(r1*rho23);

    lconn[2][0][0] = -a2*r1*s2th*irho23_dthdx2;
    lconn[2][0][1] = r1*lconn[2][0][0];
    lconn[2][0][2] = 0.;
    lconn[2][0][3] = a*r1*(a2+r2)*s2th*irho23_dthdx2;

    lconn[2][1][0] = lconn[2][0][1];
    lconn[2][1][1] = r2*lconn[2][0][0];
    lconn[2][1][2] = r2*irho2;
    lconn[2][1][3] = (a*r1*cth*sth*(r3*(2. + r1) + a2*(2.*r1*(1. + r1)*cth2 + a2*cth4 + 2*r1sth2)))*irho23_dthdx2;
    
    lconn[2][2][0] = lconn[2][0][2];
    lconn[2][2][1] = lconn[2][1][2];
    lconn[2][2][2] = -a2*cth*sth*dthdx2*irho2 + d2thdx22/dthdx2;
    lconn[2][2][3] = 0.;
    
    lconn[2][3][0] = lconn[2][0][3];
    lconn[2][3][1] = lconn[2][1][3];
    lconn[2][3][2] = lconn[2][2][3];
    lconn[2][3][3] = -cth*sth*(rho23 + a2sth2*rho2*(r1*(4. + r1) + a2cth2) + 2.*r1*a4*sth4)*irho23_dthdx2;
    
    lconn[3][0][0] = a*fac1_rho23;
    lconn[3][0][1] = r1*lconn[3][0][0];
    lconn[3][0][2] = -2.*a*r1*cth*dthdx2/(sth*rho22);
    lconn[3][0][3] = -a2sth2*fac1_rho23;

    lconn[3][1][0] = lconn[3][0][1];
    lconn[3][1][1] = a*r2*fac1_rho23;
    lconn[3][1][2] = -2*a*r1*(a2 + 2*r1*(2. + r1) + a2*c2th)*cth*dthdx2/(sth*fac2*fac2);
    lconn[3][1][3] = r1*(r1*rho22 - a2sth2*fac1)*irho23;
    
    lconn[3][2][0] = lconn[3][0][2];
    lconn[3][2][1] = lconn[3][1][2];
    lconn[3][2][2] = -a*r1*dthdx22*irho2;
    lconn[3][2][3] = dthdx2*(0.25*fac2*fac2*cth/sth + a2*r1*s2th)*irho22;
    
    lconn[3][3][0] = lconn[3][0][3];
    lconn[3][3][1] = lconn[3][1][3];
    lconn[3][3][2] = lconn[3][2][3];
    lconn[3][3][3] = (-a*r1sth2*rho22 + a3*sth4*fac1)*irho23;
    
	
}




/*
	HARM model specification routines 
*/


void interp_fourv(double X[NDIM], double ****fourv, double Fourv[NDIM]) ;
double interp_scalar(double X[NDIM], double ***var) ;

void init_model(char *args[])
{
	void init_harm3d_grid(char *);
	void init_harm3d_data(char *);

	fprintf(stderr, "reading data header...\n");
	/* Read in header and allocate space for grid data */
	init_harm3d_grid(args[3]);
	fprintf(stderr, "success\n");

	/* find dimensional quantities from black hole
		mass and its accretion rate */
	set_units(args[4]);

	fprintf(stderr, "reading data...\n");
	/* Read in the grid data */
	init_harm3d_data(args[3]);
	fprintf(stderr, "success\n");

	/* pre-compute densities, field strengths, etc. */
	init_physical_quantities() ;

	/* horizon radius */
	Rh = 1 + sqrt(1. - a * a) ;

}

/* 

	these supply basic model data to grmonty

*/

void get_model_ucov(double X[NDIM], double Ucov[NDIM])
{
	double gcov[NDIM][NDIM];

	gcov_func(X, gcov);

	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   
	   	/* sensible default value */
		Ucov[0] = -1./sqrt(-gcov[0][0]) ;
		Ucov[1] = 0. ;
		Ucov[2] = 0. ;
		Ucov[3] = 0. ;

		return ;
	}

	interp_fourv(X, ucov, Ucov) ;

}

void get_model_ucon(double X[NDIM], double Ucon[NDIM])
{

	double gcov[NDIM][NDIM] ;
	double gcon[NDIM][NDIM] ;
	double tmp[NDIM] ;

	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   	/* sensible default value */
	   	gcov_func(X, gcov) ;

		tmp[0] = -1./sqrt(-gcov[0][0]) ;
		tmp[1] = 0. ;
		tmp[2] = 0. ;
		tmp[3] = 0. ;

	   	gcon_func(gcov, gcon) ;
		Ucon[0] = 
			tmp[0]*gcon[0][0] +
			tmp[1]*gcon[0][1] +
			tmp[2]*gcon[0][2] +
			tmp[3]*gcon[0][3] ;
		Ucon[1] = 
			tmp[0]*gcon[1][0] +
			tmp[1]*gcon[1][1] +
			tmp[2]*gcon[1][2] +
			tmp[3]*gcon[1][3] ;
		Ucon[2] = 
			tmp[0]*gcon[2][0] +
			tmp[1]*gcon[2][1] +
			tmp[2]*gcon[2][2] +
			tmp[3]*gcon[2][3] ;
		Ucon[3] = 
			tmp[0]*gcon[3][0] +
			tmp[1]*gcon[3][1] +
			tmp[2]*gcon[3][2] +
			tmp[3]*gcon[3][3] ;
	
		return ;
	}
	   
	interp_fourv(X, ucon, Ucon) ;
}

void get_model_bcov(double X[NDIM], double Bcov[NDIM])
{
	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {

	   	Bcov[0] = 0. ;
	   	Bcov[1] = 0. ;
	   	Bcov[2] = 0. ;
	   	Bcov[3] = 0. ;

		return ;
	}
	interp_fourv(X, bcov, Bcov) ;
}

void get_model_bcon(double X[NDIM], double Bcon[NDIM])
{
	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {

	   	Bcon[0] = 0. ;
	   	Bcon[1] = 0. ;
	   	Bcon[2] = 0. ;
	   	Bcon[3] = 0. ;

		return ;
	}
	interp_fourv(X, bcon, Bcon) ;
}

double get_model_thetae(double X[NDIM])
{
	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   	return(0.) ;
	}
	
	return(interp_scalar(X, thetae)) ;
}


double get_model_uu(double X[NDIM])
{
	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   	return(0.) ;
	}
	
	return(interp_scalar(X, uu)) ;
}


//b field strength in Gauss
double get_model_b(double X[NDIM])
{

	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   	return(0.) ;
	}

	return(interp_scalar(X, b)) ;
}

double get_model_ne(double X[NDIM])
{
	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   	return(0.) ;
	}
	
	return(interp_scalar(X, ne)) ;
}


/** HARM utilities **/

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]) ;

/********************************************************************

				Interpolation routines

 ********************************************************************/

/* return fluid four-vector in simulation units */
void interp_fourv(double X[NDIM], double ****fourv, double Fourv[NDIM]){
	double del[NDIM],b1,b2,b3,d1,d2,d3,d4;
	int i, j, k, ip1, jp1, kp1;

	/* find the current zone location and offsets del[0], del[1] */
	Xtoijk(X, &i, &j, &k, del);

	ip1 = i + 1;
	jp1 = j + 1;
	kp1 = k + 1;
	
	//conditions at the x3 periodic boundary (at the last active zone)
	if(k==(N3-1)) kp1=0;   

	b1 = 1.-del[1];
	b2 = 1.-del[2];
	b3 = 1.-del[3];

	d1 = b1*b2;
	d3 = del[1] * b2;
	d2 = del[2] * b1;
	d4 = del[1] * del[2];

	/* Interpolate along x1,x2 first */
	Fourv[0] = d1*fourv[i][j][k][0] + d2*fourv[i][jp1][k][0] + d3*fourv[ip1][j][k][0] + d4*fourv[ip1][jp1][k][0];
	Fourv[1] = d1*fourv[i][j][k][1] + d2*fourv[i][jp1][k][1] + d3*fourv[ip1][j][k][1] + d4*fourv[ip1][jp1][k][1];
	Fourv[2] = d1*fourv[i][j][k][2] + d2*fourv[i][jp1][k][2] + d3*fourv[ip1][j][k][2] + d4*fourv[ip1][jp1][k][2];
	Fourv[3] = d1*fourv[i][j][k][3] + d2*fourv[i][jp1][k][3] + d3*fourv[ip1][j][k][3] + d4*fourv[ip1][jp1][k][3];

	/* Now interpolate above in x3 */
	Fourv[0] = b3*Fourv[0] + del[3]*(d1*fourv[i][j][kp1][0] + d2*fourv[i][jp1][kp1][0] + d3*fourv[ip1][j][kp1][0] + d4*fourv[ip1][jp1][kp1][0]);
	Fourv[1] = b3*Fourv[1] + del[3]*(d1*fourv[i][j][kp1][1] + d2*fourv[i][jp1][kp1][1] + d3*fourv[ip1][j][kp1][1] + d4*fourv[ip1][jp1][kp1][1]);
	Fourv[2] = b3*Fourv[2] + del[3]*(d1*fourv[i][j][kp1][2] + d2*fourv[i][jp1][kp1][2] + d3*fourv[ip1][j][kp1][2] + d4*fourv[ip1][jp1][kp1][2]);
	Fourv[3] = b3*Fourv[3] + del[3]*(d1*fourv[i][j][kp1][3] + d2*fourv[i][jp1][kp1][3] + d3*fourv[ip1][j][kp1][3] + d4*fourv[ip1][jp1][kp1][3]);

}

/* return	 scalar in cgs units */
double interp_scalar(double X[NDIM], double ***var)
{
	double del[NDIM],b1,b2,interp;
	int i, j, k, ip1, jp1, kp1;

	/* find the current zone location and offsets del[0], del[1] */
	Xtoijk(X, &i, &j, &k, del);

	ip1 = i+1;
	jp1 = j+1;
	kp1 = k+1;
	if(k==(N3-1)) kp1=0;   	

	b1 = 1.-del[1];
	b2 = 1.-del[2];

	/* Interpolate in x1,x2 first */

	interp = var[i][j][k]*b1*b2 + 
	  var[i][jp1][k]*b1*del[2] + 
	  var[ip1][j][k]*del[1]*b2 + 
	  var[ip1][jp1][k]*del[1]*del[2];


	/* Now interpolate above in x3 */

	interp = (1.-del[3])*interp + 
     	  del[3]*(var[i  ][j  ][kp1]*b1*b2 +
		  var[i  ][jp1][kp1]*del[2]*b1 +
		  var[ip1][j  ][kp1]*del[1]*b2 +
		  var[ip1][jp1][kp1]*del[1]*del[2]);
	

	return(interp);

}

/***********************************************************************************

					End interpolation routines

 ***********************************************************************************/

//check if this include the cuts
void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM])
{
  double phi;

	/* Map X[3] into sim range, assume startx[3] = 0 */
	phi = fmod(X[3], stopx[3]);
	//fold it to be positive and find index
	if(phi < 0.0) phi = stopx[3]+phi;
	
	//give index of a zone - zone index is moved to the grid zone center/
	*i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
	*j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
	*k = (int) ((phi  - startx[3]) / dx[3] - 0.5 + 1000) - 1000;	


	if(*i < 0) {
	  *i = 0 ;
	  del[1] = 0. ;
	}
	else if(*i > N1-2) {
	  *i = N1-2 ;
	  del[1] = 1. ;
	}
	else {
	  del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
	}
	
	if(*j < 0) {
          *j = 0 ;
          del[2] = 0. ;
        }
        else if(*j > N2-2) {
          *j = N2-2 ;
          del[2] = 1. ;
        }
        else {
          del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
        }

        if(*k < 0) {
          *k = 0 ;
          del[3] = 0. ;
        }
        else if(*k > N3-1) {
          *k = N3-1;
          del[3] = 1. ;
        }
        else {
          del[3] = (phi - ((*k + 0.5) * dx[3] + startx[3])) / dx[3];
        }

	return;
}

/* return boyer-lindquist coordinate of point */
void bl_coord(double *X, double *r, double *th)
{
	*r = exp(X[1]) + R0 ;
	*th = M_PI * X[2] + hslope*sin(2. * M_PI * X[2]);
	return;
}

void coord(int i, int j, int k, double *X)
{

	/* returns zone-centered values for coordinates */
	X[0] = startx[0];
	X[1] = startx[1] + (i + 0.5) * dx[1];
	X[2] = startx[2] + (j + 0.5) * dx[2];
	X[3] = startx[3] + (k + 0.5) * dx[3];

	return;
}


//this is set in decs.h
void set_units(char *munitstr)
{
	double MBH ;

#if(SOURCE_SGRA)
    MBH = MSGRA;
#endif
#if(SOURCE_M87)
    MBH = MM87;
#endif
#if(SOURCE_DABHB)
    MBH = MABHB;
#endif
#if(SOURCE_NT)
    MBH = 10;
#endif

	sscanf(munitstr,"%lf",&M_unit) ;

	MBH *= MSUN ;

	/** from this, calculate units of length, time, mass,
	    and derivative units **/
	L_unit = GNEWT * MBH / (CL * CL);
	T_unit = L_unit / CL;
	RHO_unit = M_unit / pow(L_unit, 3);
	U_unit = RHO_unit * CL * CL;
	B_unit = CL * sqrt(4.*M_PI*RHO_unit);

	fprintf(stderr,"L,T,M units: %g [cm] %g [sec] %g [g]\n",L_unit,T_unit,M_unit) ;
	fprintf(stderr,"rho,u,B units: %g [g cm^-3] %g [g cm^-1 sec^-2] %g [G] \n",RHO_unit,U_unit,B_unit) ;
}



void init_physical_quantities(void)
{
	int i, j, k;
        double bsq,Thetae_unit,sigma_m,beta,b2,trat;
	//	double game,gamp,ue,Te,Tp;
	double r,th,X[NDIM];
	
	for (i = 0; i < N1; i++) {
	  X[1] = startx[1] + ( i + 0.5)*dx[1];
	  for (j = 0; j < N2; j++) {
	    X[2] = startx[2] + (j+0.5)*dx[2];
	    bl_coord(X, &r, &th);

	    for (k = 0; k < N3; k++) {
	      
	      uu[i][j][k] = p[UU][i][j][k];

	      bsq= bcon[i][j][k][0] * bcov[i][j][k][0] +
		bcon[i][j][k][1] * bcov[i][j][k][1] +
		bcon[i][j][k][2] * bcov[i][j][k][2] +
		bcon[i][j][k][3] * bcov[i][j][k][3] ;
	      b[i][j][k] = sqrt(bsq)*B_unit ;
	      sigma_m=bsq/p[KRHO][i][j][k] ;

	      /* beta presciption Moscibrodzka et al. 2016, 2017*/
	      beta=p[UU][i][j][k]*(gam-1.)/0.5/bsq;
	      b2=pow(beta,2);
	      trat = trat_d * b2/(1. + b2) + trat_j /(1. + b2);
	      //Thetae_unit= (MP / ME) * 2./(2.+trat) * (gam-1.);
	      Thetae_unit = (gam - 1.) * (MP / ME) / trat;
	      thetae[i][j][k] = (p[UU][i][j][k]/p[KRHO][i][j][k])* Thetae_unit;

	      /* floor and ceiling on temp */
	      thetae[i][j][k]=MAX(thetae[i][j][k],THETAE_MIN);
	      thetae[i][j][k]=MIN(thetae[i][j][k],THETAE_MAX);

	      ne[i][j][k] = p[KRHO][i][j][k] * RHO_unit/(MP+ME) ;
	      if(thetae[i][j][k] > 100.0) ne[i][j][k]=0.0;	     
	      
	      //apply sigma cut
	      if(sigma_m > sigma_cut) {
		  ne[i][j][k]=0.0;
		  thetae[i][j][k]=0.0;
		  b[i][j][k]=0.0;
	      }
	      
	      
	    }
	  }
	}




	return ;
}


void *malloc_rank1(int n1, int size)
{
	void *A;

	if((A = malloc(n1*size)) == NULL){
		fprintf(stderr,"malloc failure in malloc_rank1\n");
		exit(123);
	}

	return A;
}

double **malloc_rank2(int n1, int n2)
{

	double **A;
	double *space;
	int i;

	space = malloc_rank1(n1*n2, sizeof(double));

	A = malloc_rank1(n1, sizeof(double *));

	for(i = 0; i < n1; i++) A[i] = &(space[i*n2]);

	return A;
}


double ***malloc_rank3(int n1, int n2, int n3)
{

	double ***A;
	double *space;
	int i,j;

	space = malloc_rank1(n1*n2*n3, sizeof(double));

	A = malloc_rank1(n1, sizeof(double *));

	for(i = 0; i < n1; i++){
		A[i] = malloc_rank1(n2,sizeof(double *));
		for(j = 0; j < n2; j++){
			A[i][j] = &(space[n3*(j + n2*i)]);
		}
	}

	return A;
}


double ****malloc_rank4(int n1, int n2, int n3, int n4)
{

	double ****A;
	double *space;
	int i,j,k;

	space = malloc_rank1(n1*n2*n3*n4, sizeof(double));

	A = malloc_rank1(n1, sizeof(double *));
	
	for(i=0;i<n1;i++){
		A[i] = malloc_rank1(n2,sizeof(double *));
		for(j=0;j<n2;j++){
			A[i][j] = malloc_rank1(n3,sizeof(double *));
			for(k=0;k<n3;k++){
				A[i][j][k] = &(space[n4*(k + n3*(j + n2*i))]);
			}
		}
	}

	return A;
}

double *****malloc_rank5(int n1, int n2, int n3, int n4, int n5)
{

	double *****A;
	double *space;
	int i,j,k,l;

	space = malloc_rank1(n1*n2*n3*n4*n5, sizeof(double));

	A = malloc_rank1(n1, sizeof(double *));
	
	for(i=0;i<n1;i++){
		A[i] = malloc_rank1(n2, sizeof(double *));
		for(j=0;j<n2;j++){
			A[i][j] = malloc_rank1(n3, sizeof(double *));
			for(k=0;k<n3;k++){
				A[i][j][k] = malloc_rank1(n4, sizeof(double *));
				for(l=0;l<n4;l++){
					A[i][j][k][l] = &(space[n5*(l + n4*(k + n3*(j + n2*i)))]);
				}
			}
		}
	}

	return A;
}

void init_storage(void)
{
	int i;

	bcon = malloc_rank4(N1,N2,N3,NDIM);
	bcov = malloc_rank4(N1,N2,N3,NDIM);
	ucon = malloc_rank4(N1,N2,N3,NDIM);
	ucov = malloc_rank4(N1,N2,N3,NDIM);
	p = (double ****)malloc_rank1(NPRIM,sizeof(double *));
	for(i = 0; i < NPRIM; i++) p[i] = malloc_rank3(N1,N2,N3);
	ne = malloc_rank3(N1,N2,N3);
	uu = malloc_rank3(N1,N2,N3);
	thetae = malloc_rank3(N1,N2,N3);
	b = malloc_rank3(N1,N2,N3);

	return;
}

/* HDF5 v1.6 API */
//#include <H5LT.h>
/* HDF5 v1.8 API */
#include <hdf5.h>
#include <hdf5_hl.h>

void init_harm3d_grid(char *fname)
{
	hid_t file_id;
	double th_end,th_cutout,t;

	file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(file_id < 0){
		fprintf(stderr,"file %s does not exist, aborting...\n",fname);
		exit(1234);
	}

	H5LTread_dataset_double(file_id,		"/Header/Grid/t", 	&t);
	H5LTread_dataset_int(file_id,		"/Header/Grid/totalsize1", 	&N1);
	H5LTread_dataset_int(file_id,		"/Header/Grid/totalsize2", 	&N2);
	H5LTread_dataset_int(file_id,		"/Header/Grid/totalsize3", 	&N3);
	H5LTread_dataset_double(file_id, 	"/Header/Grid/dx1", 		&(dx[1]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/dx2", 		&(dx[2]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/dx3", 		&(dx[3]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/startx0", 	&(startx[0]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/startx1", 	&(startx[1]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/startx2", 	&(startx[2]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/startx3", 	&(startx[3]));
	H5LTread_dataset_double(file_id, 	"/Header/Grid/a",		&a);
	H5LTread_dataset_double(file_id,	"/Header/Grid/gam",		&gam);
	H5LTread_dataset_double(file_id,	"/Header/Grid/R0",		&R0);
	H5LTread_dataset_double(file_id,	"/Header/Grid/Rin",		&Rin);
	H5LTread_dataset_double(file_id, 	"/Header/Grid/Rout",		&Rout);
	H5LTread_dataset_double(file_id, 	"/Header/Grid/h_slope",		&hslope);
	H5LTread_dataset_double(file_id,	"/Header/Grid/th_cutout",  	&th_cutout);

	/* Account for 3 ghost zones */
	startx[1] += 3*dx[1];
	startx[2] += 3*dx[2];
        startx[3] += 3*dx[3];

	stopx[0] = 1.;
	stopx[1] = startx[1]+N1*dx[1];
	stopx[2] = startx[2]+N2*dx[2];
	stopx[3] = startx[3]+N3*dx[3];

	th_beg=th_cutout;
	th_end=M_PI-th_cutout;
	th_len = th_end-th_beg;
	
	fprintf(stdout,"size: %d %d %d t= %g M\n",N1,N2,N3, t);
	fprintf(stdout,"---> a,R0 %g %g \n",a,R0);
	fprintf(stdout,"start: %g %g %g \n",startx[1],startx[2],startx[3]);
	fprintf(stdout,"th_cutout: %g  %d x %d x %d\n",th_cutout,N1,N2,N3);
	fprintf(stdout,"stop: %g %g %g \n",stopx[1],stopx[2],stopx[3]);

	init_storage();
	H5Fclose(file_id);
}

void init_harm3d_data(char *fname)
{
	hid_t file_id;
	int i,j,k,l,m;
	double X[NDIM],UdotU,ufac,udotB;
	double gcov[NDIM][NDIM],gcon[NDIM][NDIM], g;
	double dMact, Ladv, MBH;
	double r,th;
	FILE *fp;
	double BSQ,Trphi,dotJ;
	int dotJ_count;


#if(SOURCE_SGRA)
    MBH = MSGRA;
#endif
#if(SOURCE_M87)
    MBH = MM87;
#endif
#if(SOURCE_DABHB)
    MBH = MABHB;
#endif
#if(SOURCE_NT)
    MBH = 10;
#endif

	MBH=MBH*MSUN;
        L_unit = GNEWT * MBH / (CL * CL);
        T_unit = L_unit / CL;
	
	file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(file_id < 0){
		fprintf(stderr,"file %s does not exist, aborting...\n",fname);
		exit(12345);
	}

	H5LTread_dataset_double(file_id, "rho",	&(p[KRHO][0][0][0]));
	H5LTread_dataset_double(file_id, "uu",	&(p[UU][0][0][0]));
	H5LTread_dataset_double(file_id, "v1",	&(p[U1][0][0][0]));
	H5LTread_dataset_double(file_id, "v2",	&(p[U2][0][0][0]));
	H5LTread_dataset_double(file_id, "v3",	&(p[U3][0][0][0]));
	H5LTread_dataset_double(file_id, "B1",	&(p[B1][0][0][0]));
	H5LTread_dataset_double(file_id, "B2",	&(p[B2][0][0][0]));
	H5LTread_dataset_double(file_id, "B3",	&(p[B3][0][0][0]));

	H5Fclose(file_id);
	X[0] = 0.;
	X[3] = 0.;

	//fprintf(stderr,"reconstructing 4-vectors...\n");
	dMact = Ladv = 0.;
	dotJ=0.0;
	
	//reconstruction of variables at the zone center!
	for(i = 0; i < N1; i++){
	  X[1] = startx[1] + ( i + 0.5)*dx[1];
	  for(j = 0; j < N2; j++){
	    X[2] = startx[2] + (j+0.5)*dx[2];
	    gcov_func(X, gcov); // in system with cut off
	    gcon_func(gcov, gcon);
	    g = gdet_func(gcov);
	    bl_coord(X, &r, &th);

	    for(k = 0; k < N3; k++){
	      UdotU = 0.;
	      X[3] = startx[3] + (k+0.5)*dx[3];
	      
	      //the four-vector reconstruction should have gcov and gcon and gdet using the modified coordinates
	      //interpolating the four vectors to the zone center !!!!
	      for(l = 1; l < NDIM; l++) for(m = 1; m < NDIM; m++) UdotU += gcov[l][m]*p[U1+l-1][i][j][k]*p[U1+m-1][i][j][k];
	      ufac = sqrt(-1./gcon[0][0]*(1 + fabs(UdotU)));
	      ucon[i][j][k][0] = -ufac*gcon[0][0];
	      for(l = 1; l < NDIM; l++) ucon[i][j][k][l] = p[U1+l-1][i][j][k] - ufac*gcon[0][l];
	      lower(ucon[i][j][k], gcov, ucov[i][j][k]);

	      //reconstruct the magnetic field three vectors
	      udotB = 0.;
	      for(l = 1; l < NDIM; l++) udotB += ucov[i][j][k][l]*p[B1+l-1][i][j][k];
	      bcon[i][j][k][0] = udotB;
	      for(l = 1; l < NDIM; l++) bcon[i][j][k][l] = (p[B1+l-1][i][j][k] + ucon[i][j][k][l]*udotB)/ucon[i][j][k][0];
	      lower(bcon[i][j][k], gcov, bcov[i][j][k]);

	      if(i <= 20) dMact += g * p[KRHO][i][j][k] * ucon[i][j][k][1] ;
	      if(i >= 20 && i < 40) Ladv += g * p[UU][i][j][k] * ucon[i][j][k][1] * ucov[i][j][k][0] ;
	      
	    }
	  }
	}

	dMact *= dx[3]*dx[2] ;
	dMact /= 21. ;
	Ladv *= dx[3]*dx[2] ;
	Ladv /= 21. ;

	fprintf(stderr,"dMact: %g [code]\n",dMact) ;
	fprintf(stderr,"Ladv: %g [code]\n",Ladv) ;
	fprintf(stderr,"Mdot: %g [g/s] \n",-dMact*M_unit/T_unit) ;
	fprintf(stderr,"Mdot: %g [MSUN/YR] \n",-dMact*M_unit/T_unit/(MSUN / YEAR)) ;
        double Mdotedd=4.*M_PI*GNEWT*MBH*MP/CL/0.1/SIGMA_THOMSON;
        fprintf(stderr,"Mdot: %g [Mdotedd]\n",-dMact*M_unit/T_unit/Mdotedd) ;
        fprintf(stderr,"Mdotedd: %g [g/s]\n",Mdotedd) ;
        fprintf(stderr,"Mdotedd: %g [MSUN/YR]\n",Mdotedd/(MSUN/YEAR)) ;

}

/* this is from illinois version that does not need second derivative of theta, works for funky cooridnates or whatever */
double theta_func(double x[NDIM])
{
    double r,th;
    bl_coord(x, &r, &th);
    return th;
}

double root_find2(double X[NDIM])
{
    double th = X[2];
    double thb, thc;

    double Xa[NDIM], Xb[NDIM], Xc[NDIM];
    Xa[1] = log(X[1]);
    Xa[3] = X[3];
    Xb[1] = Xa[1];
    Xb[3] = Xa[3];
    Xc[1] = Xa[1];
    Xc[3] = Xa[3];

    if (X[2] < M_PI / 2.) {
        Xa[2] = 0. - SMALL;
        Xb[2] = 0.5 + SMALL;
    } else {
        Xa[2] = 0.5 - SMALL;
        Xb[2] = 1. + SMALL;
    }

    thb = theta_func(Xb);

    /* bisect for a bit */
    double tol = 1.e-9;
    for (int i = 0; i < 1000; i++) {
        Xc[2] = 0.5 * (Xa[2] + Xb[2]);
        thc = theta_func(Xc);

        if ((thc - th) * (thb - th) < 0.)
            Xa[2] = Xc[2];
        else
            Xb[2] = Xc[2];

        double err = theta_func(Xc) - th;
        if (fabs(err) < tol)
            break;
    }

    return (Xa[2]);

}

