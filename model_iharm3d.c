#include "decs.h"
#include <string.h>

/*

  first geometry used in the model and then 
	HARM model specification routines 
 */

/* geometry used in this model*/


/*


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
//this seems correct
#define MUNULOOP for(int mu=0;mu<NDIM;mu++) for(int nu=0;nu<NDIM;nu++)
void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{
    // Jacobian with respect to KS basis where X is given in
    // non-KS basis
  MUNULOOP dxdX[mu][nu] = 0.;

  dxdX[0][0] = 1.;
  dxdX[1][1] = exp(X[1]);
  dxdX[2][1] = -exp(mks_smooth * (startx[1] - X[1])) * mks_smooth
      * (
	  M_PI / 2. -
	  M_PI * X[2]
	  + poly_norm * (2. * X[2] - 1.)
	  * (1
	     + (pow((-1. + 2 * X[2]) / poly_xt, poly_alpha))
	     / (1 + poly_alpha))
	  - 1. / 2. * (1. - hslope) * sin(2. * M_PI * X[2]));
  dxdX[2][2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2])
      + exp(mks_smooth * (startx[1] - X[1]))
      * (-M_PI
	 + 2. * poly_norm
	 * (1.
	    + pow((2. * X[2] - 1.) / poly_xt, poly_alpha)
	    / (poly_alpha + 1.))
	 + (2. * poly_alpha * poly_norm * (2. * X[2] - 1.)
	    * pow((2. * X[2] - 1.) / poly_xt, poly_alpha - 1.))
	 / ((1. + poly_alpha) * poly_xt)
	 - (1. - hslope) * M_PI * cos(2. * M_PI * X[2]));
  dxdX[3][3] = 1.;
  
  
}


/* 
   Current metric: modified Kerr-Schild, squashed in theta
   to give higher resolution at the equator 
*/
    
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

void get_connection(double X[NDIM], double lconn[NDIM][NDIM][NDIM])
{

    //for these coordinates it is safer to go numerically
    get_connection_num(X,lconn);

}


/*
double ****bcon;
double ****bcov;
double ****ucon;
double ****ucov;
double ****p;
double ***ne;
double ***uu;
double ***thetae;
double ***b;
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

	//get_model_ucon(X, Ucon);
	//lower(Ucon, gcov, Ucov);

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

	 // since we read from data, adjust i,j,k for ghost zones
	i += 1;
	j += 1;
	k += 1;
	
	ip1 = i + 1;
	jp1 = j + 1;
	kp1 = k + 1;
	
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
	//new

	//no interpolation of vectors at all
	/*
	Fourv[0]=fourv[i][j][k][0];
	Fourv[1]=fourv[i][j][k][1];
	Fourv[2]=fourv[i][j][k][2];
	Fourv[3]=fourv[i][j][k][3];
	*/
}

/* return	 scalar in cgs units */
double interp_scalar(double X[NDIM], double ***var)
{
	double del[NDIM],b1,b2,interp;
	int i, j, k, ip1, jp1, kp1;

	/* find the current zone location and offsets del[0], del[1] */
	Xtoijk(X, &i, &j, &k, del);

	// since we read from data, adjust i,j,k for ghost zones
	i += 1;
	j += 1;
	k += 1;

	ip1 = i+1;
	jp1 = j+1;
	kp1 = k+1;

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
	
	//new, no interpolations what so ever
	//interp=var[i][j][k];
	/* use bilinear interpolation to find rho; piecewise constant
	   near the boundaries */
	
	return(interp);

}

/***********************************************************************************

					End interpolation routines

 ***********************************************************************************/


void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM])
{
  double phi;

  phi = fmod(X[3], stopx[3]);
  if(phi < 0.0) phi = stopx[3]+phi;

  // get provisional zone index. see note above function for details. note we
  // shift to zone centers because that's where variables are most exact.
  *i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
  *j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
  *k = (int) ((phi  - startx[3]) / dx[3] - 0.5 + 1000) - 1000;

  // exotic coordinate systems sometime have issues. use this block to enforce
  // reasonable limits on *i,*j and *k. in the normal coordinate systems, this
  // block should never fire.
  if (*i < -1) *i = -1;
  if (*j < -1) *j = -1;
  if (*k < -1) *k = -1;
  if (*i >= N1) *i = N1-1;
  if (*j >= N2) *j = N2-1;
  if (*k >= N3) *k = N3-1;

  // now construct del
  del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
  del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
  del[3] = (phi - ((*k + 0.5) * dx[3] + startx[3])) / dx[3];

  // and enforce limits on del (for exotic coordinate systems)
  for (int i=0; i<4; ++i) {
    if (del[i] < 0.) del[i] = 0.;
    if (del[i] >= 1.) del[i] = 1.;
  }
	  
	return;
}

//#define SINGSMALL (1.E-20)
/* return boyer-lindquist coordinate of point */
void bl_coord(double *X, double *r, double *th)
{

     *r = exp(X[1])  ;

     double thG = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
     double y = 2 * X[2] - 1.;
     double thJ = poly_norm * y
	 * (1. + pow(y / poly_xt, poly_alpha) / (poly_alpha + 1.)) + 0.5 * M_PI;
     *th = thG + exp(mks_smooth * (startx[1] - X[1])) * (thJ - thG);
        
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
    MBH = 1;
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
        double bsq,Thetae_unit,beta,b2,trat,sigma_m;

	//all zones: active+ghost
	for (i = 0; i < N1+2; i++) {
	  for (j = 0; j < N2+2; j++) {
	    for (k = 0; k < N3+2; k++) {
	      
	      uu[i][j][k] = p[UU][i][j][k];
	      bsq= bcon[i][j][k][0] * bcov[i][j][k][0] +
		bcon[i][j][k][1] * bcov[i][j][k][1] +
		bcon[i][j][k][2] * bcov[i][j][k][2] +
		bcon[i][j][k][3] * bcov[i][j][k][3] ;
	      b[i][j][k] = sqrt(bsq)*B_unit ;
	      sigma_m=bsq/p[KRHO][i][j][k] ;
	      beta=p[UU][i][j][k]*(gam-1.)/0.5/bsq;
	      b2=beta*beta;
	      trat = trat_d * b2/(1. + b2) + trat_j /(1. + b2);
	      Thetae_unit = (MP/ME) * (game-1.) * (gamp-1.) / ( (gamp-1.) + (game-1.)*trat );
	      thetae[i][j][k] = (p[UU][i][j][k]/p[KRHO][i][j][k])* Thetae_unit;
	      thetae[i][j][k] = MAX(thetae[i][j][k], THETAE_MIN);
	      ne[i][j][k] = p[KRHO][i][j][k] * RHO_unit/(MP+ME) ;
	      	      
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

	bcon = malloc_rank4(N1+2,N2+2,N3+2,NDIM);
	bcov = malloc_rank4(N1+2,N2+2,N3+2,NDIM);
	ucon = malloc_rank4(N1+2,N2+2,N3+2,NDIM);
	ucov = malloc_rank4(N1+2,N2+2,N3+2,NDIM);
	p = (double ****)malloc_rank1(NPRIM,sizeof(double *));
	for(i = 0; i < NPRIM; i++) p[i] = malloc_rank3(N1+2,N2+2,N3+2);
	ne = malloc_rank3(N1+2,N2+2,N3+2);
	uu = malloc_rank3(N1+2,N2+2,N3+2);
	thetae = malloc_rank3(N1+2,N2+2,N3+2);
	b = malloc_rank3(N1+2,N2+2,N3+2);

	return;
}

/* HDF5 v1.6 API */
//#include <H5LT.h>

/* HDF5 v1.8 API */
#include <hdf5.h>
#include <hdf5_hl.h>

/* Harm3d globals */
/*
extern double ****bcon;
extern double ****bcov;
extern double ****ucon;
extern double ****ucov;
extern double ****p;
extern double ***ne;
extern double ***uu;
extern double ***thetae;
extern double ***b;
*/

void init_harm3d_grid(char *fname)
{
	hid_t file_id;
	double th_cutout,t;

	file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(file_id < 0){
		fprintf(stderr,"file %s does not exist, aborting...\n",fname);
		exit(1234);
	}

	H5LTread_dataset_int(file_id,		"/header/n1", 	&N1);
	H5LTread_dataset_int(file_id,		"/header/n2", 	&N2);
	H5LTread_dataset_int(file_id,		"/header/n3", 	&N3);
	H5LTread_dataset_double(file_id,	"/header/gam",		&gam);
	H5LTread_dataset_double(file_id,	"/header/gam_e",		&game);
	H5LTread_dataset_double(file_id,	"/header/gam_p",		&gamp);		
	H5LTread_dataset_double(file_id, 	"/header/geom/startx1", 	&(startx[1]));
	H5LTread_dataset_double(file_id, 	"/header/geom/startx2", 	&(startx[2]));
	H5LTread_dataset_double(file_id, 	"/header/geom/startx3", 	&(startx[3]));
	H5LTread_dataset_double(file_id, 	"/header/geom/dx1", 		&(dx[1]));
	H5LTread_dataset_double(file_id, 	"/header/geom/dx2", 		&(dx[2]));
	H5LTread_dataset_double(file_id, 	"/header/geom/dx3", 		&(dx[3]));

      

	H5LTread_dataset_double(file_id, 	"/header/geom/mmks/a",		&a);
	//mads
	//H5LTread_dataset_double(file_id,	"/header/geom/mmks/r_in",	&Rin);
	//H5LTread_dataset_double(file_id, 	"/header/geom/mmks/r_out",	&Rout);
	//sane
	//H5LTread_dataset_double(file_id,	"/header/geom/mmks/Rin",	&Rin);
	//H5LTread_dataset_double(file_id, 	"/header/geom/mmks/Rout",	&Rout);
	H5LTread_dataset_double(file_id, 	"/header/geom/mmks/hslope",	&hslope);	
	/*additional parameter that describe theta*/
	/* mks_smooth,poly_alpha,poly_xt */
	H5LTread_dataset_double(file_id, 	"/header/geom/mmks/mks_smooth",	&mks_smooth);
	H5LTread_dataset_double(file_id, 	"/header/geom/mmks/poly_alpha",	&poly_alpha);
	H5LTread_dataset_double(file_id, 	"/header/geom/mmks/poly_xt",	&poly_xt);	
	poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha + 1.)*1./pow(poly_xt, poly_alpha));
	

	/* for test 
	    H5LTread_dataset_double(file_id, 	"/header/geom/mks/a",		&a);	
	    H5LTread_dataset_double(file_id,	"/header/geom/mks/r_in",	&Rin);
	    H5LTread_dataset_double(file_id, 	"/header/geom/mks/r_out",	&Rout);
	    H5LTread_dataset_double(file_id, 	"/header/geom/mks/hslope",	&hslope);
	*/
	
	stopx[0] = 1.;
	stopx[1] = startx[1]+N1*dx[1];
	stopx[2] = startx[2]+N2*dx[2];
	stopx[3] = startx[3]+N3*dx[3];

	th_beg = 0.0174;
	
	fprintf(stdout,"size: %d %d %d M  gam=%g game=%g gamp=%g t=%g M mks_smooth=%g\n",N1,N2,N3,gam,game,gamp,t,mks_smooth);
	fprintf(stdout,"---> a=%g Rout=%g hslope=%g \n",a,Rout,hslope);
	fprintf(stdout,"start: %g %g %g \n",startx[1],startx[2],startx[3]);
	fprintf(stdout,"th_cutout: %g  %d x %d x %d\n",th_cutout,N1,N2,N3);
	fprintf(stdout,"stop: %g %g %g \n",stopx[1],stopx[2],stopx[3]);

	init_storage();
	H5Fclose(file_id);
}

int hdf5_read_array(void *data, const char *name, size_t rank,
		    hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims, hsize_t *mstart, hsize_t hdf5_type, char *fname)
{

    hid_t filespace = H5Screate_simple(rank, fdims, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, fstart, NULL, fcount,NULL);
    hid_t memspace = H5Screate_simple(rank, mdims, NULL);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mstart, NULL, fcount,NULL);

    hid_t file_id ;

    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

    //for test
    //file_id = H5Fopen("SANE_a+0.94_288_0900_MKS.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    
    
    hid_t dset_id = H5Dopen(file_id, "prims"); //1.6 vs 1.8
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    herr_t err = H5Dread(dset_id, hdf5_type, memspace, filespace, plist_id, data);

    H5Dclose(dset_id);
    H5Pclose(plist_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    
  return 0;
}
void init_harm3d_data(char *fname)
{

	int i,j,k,l,m;
	double X[NDIM],UdotU,ufac,udotB;
	double gcov[NDIM][NDIM],gcon[NDIM][NDIM], g;
	double dMact, Ladv, MBH;
	double r,th;
	FILE *fp;
	//double BSQ,Trphi,dotJ;
	//int dotJ_count;

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
    MBH = 1;
#endif
	

	MBH=MBH*MSUN;
        L_unit = GNEWT * MBH / (CL * CL);
        T_unit = L_unit / CL;

	hid_t file_id;
	file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(file_id < 0){
		fprintf(stderr,"file %s does not exist, aborting...\n",fname);
		exit(12345);
	}


       	fprintf(stderr,"data incoming...from --------> %s \n",fname);
	
	int n_prims;
	H5LTread_dataset_int(file_id,        "/header/n_prim",        &n_prims);
	hsize_t fdims[] = { N1, N2, N3, n_prims };
	hsize_t fstart[] = { 0, 0, 0, 0 };
	hsize_t fcount[] = { N1, N2, N3, 1 };
	hsize_t mdims[] = { N1+2, N2+2, N3+2, 1 };
	hsize_t mstart[] = { 1, 1, 1, 0 };
	 
	fstart[3] = 0;
	hdf5_read_array(p[KRHO][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE,fname);
	fstart[3] = 1;
	hdf5_read_array(p[UU][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE,fname);
	fstart[3] = 2;
	hdf5_read_array(p[U1][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE,fname);
	fstart[3] = 3;
	hdf5_read_array(p[U2][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE,fname);
	fstart[3] = 4;
	hdf5_read_array(p[U3][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE,fname);
	fstart[3] = 5;
	hdf5_read_array(p[B1][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE,fname);
	fstart[3] = 6;
	hdf5_read_array(p[B2][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE,fname);
	fstart[3] = 7;
	hdf5_read_array(p[B3][0][0], "prims", 4, fdims, fstart, fcount, mdims, mstart, H5T_IEEE_F64LE,fname); 

	H5Fclose(file_id);
	X[0] = 0.;
	X[3] = 0.;

	//fprintf(stderr,"reconstructing 4-vectors...\n");
	dMact = Ladv = 0.;
	//dotJ=0.0;
	
	//reconstruction of variables at the zone center! here ghost zone in primitive read in from the file, they are zero anyway, i think
	//in active zones
	for(i = 1; i < N1+1; i++){
	  for(j = 1; j < N2+1; j++){
	    coord(i-1,j-1,0,X); //shift for ghost zones?
	    gcov_func(X, gcov); // in system with cut off
	    gcon_func(gcov, gcon);
	    g = gdet_func(gcov);
	    bl_coord(X, &r, &th);
	    for(k = 1; k < N3+1; k++){
	      coord(i-1,j-1,k,X);
	      //the four-vector reconstruction should have gcov and gcon and gdet using the modified coordinates
	      //interpolating the four vectors to the zone center !!!!
	      UdotU = 0.;
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

	      if(i <= 21) dMact += g * p[KRHO][i][j][k] * ucon[i][j][k][1] ;
	      if(i >= 20 && i < 40 ) Ladv += g * p[UU][i][j][k] * ucon[i][j][k][1] * ucov[i][j][k][0] ;

	      //compute dimentionless T^r_phi over phi and theta
	      /*
	      if ( i == (N1-2) ){
		BSQ=bcon[i][j][k][0]*bcov[i][j][k][0]+
		  bcon[i][j][k][1]*bcov[i][j][k][1]+
		  bcon[i][j][k][2]*bcov[i][j][k][2]+
		  bcon[i][j][k][3]*bcov[i][j][k][3];
		Trphi =  (p[KRHO][i][j][k]+gam*p[UU][i][j][k]+BSQ)*ucon[i][j][k][1]*ucov[i][j][k][3]-bcon[i][j][k][1]*bcov[i][j][k][3] ;
		dotJ += Trphi*g*dx[2]*dx[3];
	      }
	      */
	      
	    }
	  }
	}


	//for r outflow
	for(j = 1; j < N2+1; ++j){
	  for(k = 1; k < N3+1; ++k){

	    for(l = 0; l < NDIM; ++l){
	      bcon[0   ][j][k][l] = bcon[1 ][j][k][l];
	      bcon[N1+1][j][k][l] = bcon[N1][j][1][l];

	      bcov[0   ][j][k][l] = bcov[1 ][j][k][l];
	      bcov[N1+1][j][k][l] = bcov[N1][j][1][l];

	      ucon[0   ][j][k][l] = ucon[1 ][j][k][l];
	      ucon[N1+1][j][k][l] = ucon[N1][j][1][l];

	      ucov[0   ][j][k][l] = ucov[1 ][j][k][l];
	      ucov[N1+1][j][k][l] = ucov[N1][j][1][l];
	    }
	    //prim
	    for(l = 0; l < NPRIM; ++l){
	      p[l][0   ][j][k]=p[l][1 ][j][k];
	      p[l][N1+1][j][k]=p[l][N1][j][k];
	    }
	    
	    
	  }
	}

	 // elevation -- flip (this is a rotation by pi)

	for (i=0; i<N1+2; ++i) {
	  for (k=1; k<N3+1; ++k) {
	    if (N3%2 == 0) {
	      int kflip = ( k + (N3/2) ) % N3;
	      for (l=0; l<NDIM; ++l) {
		bcon[i][0][k][l] = bcon[i][1][kflip][l];
		bcon[i][N2+1][k][l] = bcon[i][N2][kflip][l];
		bcov[i][0][k][l] = bcov[i][1][kflip][l];
		bcov[i][N2+1][k][l] = bcov[i][N2][kflip][l];
		ucon[i][0][k][l] = ucon[i][1][kflip][l];
		ucon[i][N2+1][k][l] = ucon[i][N2][kflip][l];
		ucov[i][0][k][l] = ucov[i][1][kflip][l];
		ucov[i][N2+1][k][l] = ucov[i][N2][kflip][l];
	      }
	      for (l=0; l<NPRIM; ++l) {
		p[l][i][0][k] = p[l][i][1][kflip];
		p[l][i][N2+1][k] = p[l][i][N2][kflip];
	      }
	    } else {
	      int kflip1 = ( k + (N3/2) ) % N3;
	      int kflip2 = ( k + (N3/2) + 1 ) % N3;
	      for (l=0; l<NDIM; ++l) {
		bcon[i][0][k][l]    = ( bcon[i][1][kflip1][l] 
					+ bcon[i][1][kflip2][l] ) / 2.;
		bcon[i][N2+1][k][l] = ( bcon[i][N2][kflip1][l]
					+ bcon[i][N2][kflip2][l] ) / 2.;
		bcov[i][0][k][l]    = ( bcov[i][1][kflip1][l]
					+ bcov[i][1][kflip2][l] ) / 2.;
		bcov[i][N2+1][k][l] = ( bcov[i][N2][kflip1][l] 
					+ bcov[i][N2][kflip2][l] ) / 2.;
		ucon[i][0][k][l]    = ( ucon[i][1][kflip1][l]
					+ ucon[i][1][kflip2][l] ) / 2.;
		ucon[i][N2+1][k][l] = ( ucon[i][N2][kflip1][l]
					+ ucon[i][N2][kflip2][l] ) / 2.;
		ucov[i][0][k][l]    = ( ucov[i][1][kflip1][l] 
					+ ucov[i][1][kflip2][l] ) / 2.;
		ucov[i][N2+1][k][l] = ( ucov[i][N2][kflip1][l] 
					+ ucov[i][N2][kflip2][l] ) / 2.;
	      }
	      for (l=0; l<NPRIM; ++l) {
		p[l][i][0][k]    = ( p[l][i][1][kflip1]
				     + p[l][i][1][kflip2] ) / 2.;
		p[l][i][N2+1][k] = ( p[l][i][N2][kflip1]
				     + p[l][i][N2][kflip2] ) / 2.;
	      }
	    }
	  }
	}
	
	//for phi, periodic
	for(i = 1; i < N1+1; i++){
	  for(j = 1; j < N2+1; ++j){

	    for(l = 0; l < NDIM; l++){
	      bcon[i][j][0   ][l] = bcon[i][j][N3][l];
	      bcon[i][j][N3+1][l] = bcon[i][j][1 ][l];

	      bcov[i][j][0   ][l] = bcov[i][j][N3][l];
	      bcov[i][j][N3+1][l] = bcov[i][j][1 ][l];

	      ucon[i][j][0   ][l] = ucon[i][j][N3][l];
	      ucon[i][j][N3+1][l] = ucon[i][j][1 ][l];

	      ucov[i][j][0   ][l] = ucov[i][j][N3][l];
	      ucov[i][j][N3+1][l] = ucov[i][j][1 ][l];
	    }
	    //prim
	    for(l = 0; l < NPRIM; ++l){
	      p[l][i][j][0   ]=p[l][i][j][N3];
	      p[l][i][j][N3+1]=p[l][i][j][1 ];
	    }
	    
	  }
	}
	

	dMact *= dx[3]*dx[2] ;
	dMact /= 21. ;
	Ladv *= dx[3]*dx[2] ;
	Ladv /= 21. ;

	/*
	double Mstar=0.4*MSUN;
	double Jorb= (Mstar*MBH)/(Mstar+MBH)*sqrt(GNEWT*(MBH+Mstar)*3.79*RSUN); //this is in cgs for sure
	fprintf(stderr,"Jorb:%g [gcm2/s]\n",Jorb) ;
	*/
	/*
	fprintf(stderr,"dotJ: %g [code]\n",dotJ) ;
	double dotJ_unit=M_unit*pow(L_unit,2)/pow(T_unit,2); 
	dotJ *= dotJ_unit; 
	fprintf(stderr,"dotJ: %g [gcm2/s2] \n",dotJ);
	exit(1);
	*/

	fprintf(stderr,"dMact: %g [code]\n",dMact) ;
	fprintf(stderr,"Ladv: %g [code]\n",Ladv) ;
	fprintf(stderr,"Mdot: %g [g/s] \n",-dMact*M_unit/T_unit) ;
	fprintf(stderr,"Mdot: %g [MSUN/YR] \n",-dMact*M_unit/T_unit/(MSUN / YEAR)) ;
        double Mdotedd=4.*M_PI*GNEWT*MBH*MP/CL/0.1/SIGMA_THOMSON;
        fprintf(stderr,"Mdot: %g [Mdotedd]\n",-dMact*M_unit/T_unit/Mdotedd) ;
        fprintf(stderr,"Mdotedd: %g [g/s]\n",Mdotedd) ;
        fprintf(stderr,"Mdotedd: %g [MSUN/YR]\n",Mdotedd/(MSUN/YEAR)) ;

}



/* check this*/
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




