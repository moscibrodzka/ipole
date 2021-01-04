#include "decs.h"


/*

	Novikov-Thorne model specification routines 

*/



/*geometrical routines */
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

double rhorizon_calc(int pos_sign)
{
  double sign;
  sign = (pos_sign) ? 1. : -1.;
  return(1. + sign*sqrt((1.-a)*(1.+a)) );
}



//from Illinois group, coefficents needed to set up metric corrections for modified coodtinates FMKS
#define MUNULOOP for(int mu=0;mu<NDIM;mu++) for(int nu=0;nu<NDIM;nu++)
void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{
    // Jacobian with respect to KS basis where X is given in
    // non-KS basis
    MUNULOOP dxdX[mu][nu] = 0.;

    //traditional coordinates with cut and with hslope
  
    dxdX[0][0] = 1.;
    dxdX[1][1] = exp(X[1]);
    dxdX[2][2] = M_PI;
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
  double gcov_bl[NDIM][NDIM];
  MUNULOOP gcov_bl[mu][nu] = 0.;

  bl_coord(X,&r,&th);

  cth = cos(th);
  sth = sin(th);
  s2 = sth * sth;
  rho2 = r * r + a * a * cth * cth;
  
  /* BL coordinates here */
  double a2 = a * a;
  double r2 = r * r;
  double DD = 1. - 2. / r + a2 / r2;
  double mu = 1. + a2 * cth * cth / r2;
  gcov_bl[0][0] = -(1. - 2. / (r * mu));
  gcov_bl[0][3] = -2. * a * s2 / (r * mu);
  gcov_bl[3][0] = gcov_bl[0][3];
  gcov_bl[1][1] = mu / DD;
  gcov_bl[2][2] = r2 * mu;
  gcov_bl[3][3] = r2 * sth * sth * (1. + a2 / r2 + 2. * a2 * s2 / (r2 * r * mu));
  
    
  // convert from ks metric to a modified one using Jacobian
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);
  MUNULOOP
    {
      gcov[mu][nu] = 0;
      for (int lam = 0; lam < NDIM; ++lam) {
	for (int kap = 0; kap < NDIM; ++kap) {
	  gcov[mu][nu] += gcov_bl[lam][kap] * dxdX[lam][mu] * dxdX[kap][nu];
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



void init_model(char *args[])
{
	fprintf(stderr, "setting units...");
	/* find dimensional quantities from black hole
		mass and its accretion rate */
	set_units(args[1]);
	Rh = 1 + sqrt(1. - a * a) ;
	Rout = 100.;

}

void bl_coord(double *X, double *r, double *th)
{   
        *r = exp(X[1])  ;
        *th = M_PI*X[2]  ;
}


void linear_interp_I_Q_disk(double mu, double *f1, double *f2);
void get_NT_S(double X[NDIM], double Kcon[NDIM], double F, double complex N_coord[NDIM][NDIM]);
double get_NT_nu(double X[NDIM], double Kcon[NDIM]);
double bnu(double nu, double T);
double krolikc(double r, double a, double r_isco);
double get_weight(double *xx, double x, int *jlo);



double bnu(double nu, double T) {
  return 2 * HPL * nu*nu*nu / (CL * CL) / (exp(HPL * nu / (KBOL * T)) - 1);
}

double krolikc(double r, double a, double r_isco)
{
  double y = sqrt(r);
  double yms = sqrt(r_isco);
  double y1 = 2. * cos(1. / 3. * (acos(a) - M_PI));
  double y2 = 2. * cos(1. / 3. * (acos(a) + M_PI));
  double y3 = -2. * cos(1. / 3. * acos(a));
  double arg1 = 3. * a / (2. * y);
  double arg2 = 3. * pow(y1 - a, 2) / (y * y1 * (y1 - y2) * (y1 - y3));
  double arg3 = 3. * pow(y2 - a, 2) / (y * y2 * (y2 - y1) * (y2 - y3));
  double arg4 = 3. * pow(y3 - a, 2) / (y * y3 * (y3 - y1) * (y3 - y2));

  return 1. - yms / y - arg1 * log(y / yms) - arg2 * log((y - y1) / (yms - y1))
      - arg3 * log((y - y2) / (yms - y2)) - arg4 * log((y - y3) / (yms - y3));
}

double get_NT_nu(double X[NDIM], double Kcon[NDIM]){

    double omega,r,th,AA,m1oA2;
    double gcov[NDIM][NDIM],gcon[NDIM][NDIM],Ucon[NDIM],Ucov[NDIM];

    bl_coord(X, &r, &th);
    gcov_func(X, gcov);
    omega=1./(pow(r,1.5)+a);
    m1oA2=gcov[0][0] + 2.*omega*gcov[0][3] + omega*omega*gcov[3][3];
    AA=sqrt(-1./m1oA2);
  
    /* Keplerian four-velocity around Kerr black hole in BL coordinates*/
    Ucon[0] = AA;
    Ucon[1] = 0.;
    Ucon[2] = 0.;
    Ucon[3] = AA*omega;
    
    //gcov and gcon for MBL coordinates
    lower(Ucon,gcov,Ucov);

    return get_fluid_nu(Kcon, Ucov);

}


void get_NT_S(double X[NDIM], double Kcon[NDIM], double F, double complex N_coord[NDIM][NDIM]){
    
    int j,k;
    double omega,r,th,AA,m1oA2;
    double gcov[NDIM][NDIM],gcon[NDIM][NDIM];
    double Ecov[NDIM][NDIM],Econ[NDIM][NDIM];
    double complex N_tetrad[NDIM][NDIM];
    double Ucon[NDIM],Ucov[NDIM];
    double SI,SQ,SU,SV,f1,f2;
    
    bl_coord(X, &r, &th);
    gcov_func(X, gcov);
    omega=1./(pow(r,1.5)+a);
    m1oA2=gcov[0][0] + 2.*omega*gcov[0][3] + omega*omega*gcov[3][3];
    AA=sqrt(-1./m1oA2);
  
    /* Keplerian four-velocity around Kerr black hole in BL coordinates, this is in both BL and MBL*/
    Ucon[0] = AA;
    Ucon[1] = 0.;
    Ucon[2] = 0.;
    Ucon[3] = AA*omega;

    //gcov and gcon for MBL coordinates
    lower(Ucon,gcov,Ucov);

    double nu=get_fluid_nu(Kcon, Ucov);
    
    //make unit vector ncov
    double ncon[NDIM],ncov[NDIM],mu;
    ncon[0]=0.0;
    ncon[1]=0.0;
    ncon[2]=1.0;
    ncon[3]=0.0;
    normalize(ncon,gcov);
    lower(ncon,gcov,ncov);
    make_plasma_tetrad(Ucon,Kcon,ncon,gcov,Econ,Ecov);
    
    mu = fabs(cos(get_bk_angle2(X, Kcon, Ucov, ncon, ncov)));
    linear_interp_I_Q_disk(mu,&f1,&f2);

//    fprintf(stdout,"%g %g \n",f1,f2);

    //f1 >> f2
    SI=F*f1; //if here f2 then QU are correct but I not, if f1 here I is fine but QU not
    SQ=F*f1*f2;
    SU=0.0;
    SV=0.0;

    SI = SI/nu/nu/nu;
    SQ = SQ/nu/nu/nu;
    
    stokes_to_tensor(SI, SQ, SU, SV, N_tetrad);
    complex_tetrad_to_coord_rank2(N_tetrad, Econ, N_coord);
    
}


double get_weight(double *xx, double x, int *jlo)
{
  //Get the value _before_ x in the table
  while (xx[*jlo] < x) {
    ++*jlo;
  }
  --*jlo;
  // Return weight for table values jlo, jlo+1
  return (x - xx[*jlo]) / (xx[*jlo + 1] - xx[*jlo]);
}

//mu is angle between k and normal to the disk plane (z axis?)
//given mu, returns I and Q
void linear_interp_I_Q_disk(double mu, double *f1, double *f2)
{
    
    double mu_tab[21]={ 0.000000,
		       0.050000,
		       0.100000,
		       0.150000,
		       0.200000,
		       0.250000,
		       0.300000,
		       0.350000,
		       0.400000,
		       0.450000,
		       0.500000,
		       0.550000,
		       0.600000,
		       0.650000,
		       0.700000,
		       0.750000,
		       0.800000,
		       0.850000,
		       0.900000,
		       0.950000,
		       1.000000};
    
    double I_tab[21]={0.414410,
		      0.474900,
		      0.523970,
		      0.570010,
		      0.614390,
		      0.657700,
		      0.700290,
		      0.742340,
		      0.783980,
		      0.825300,
		      0.866370,
		      0.907220,
		      0.947890,
		      0.988420,
		      1.02882,
		      1.06911,
		      1.10931,
		      1.14943,
		      1.18947,
		      1.22945,
		      1.26938}; 

    double Q_tab[21]={ 0.117130,
		       0.0897900,
		       0.0744800,
		       0.0631100,
		       0.0541000,
		       0.0466700,
		       0.0404100,
		       0.0350200,
		       0.0303300,
		       0.0261900,
		       0.0225200,
		       0.0192300,
		       0.0162700,
		       0.0135800,
		       0.0111230,
		       0.00888000,
		       0.00681800,
		       0.00491900,
		       0.00315500,
		       0.00152200,
		       0.00000};

    
        int i;
	double w=get_weight(mu_tab,mu,&i);
	
	//Stokes I/F
        *f1= (1. - w) * I_tab[i] + w * I_tab[i+1];
	//Stokes Q/I
	*f2= (1. - w) * Q_tab[i] + w * Q_tab[i+1];

	return ;
}



// below is used and should be here to stop integration within ipolarray.c
// if ne=0 then just transport step

/* Keplerian velocity, BL-> KS' coordinates */
void get_model_ucov(double X[NDIM], double Ucov[NDIM])
{
  int j,k;
  double omega,l,r,r2,th,r32,r12,rho2,cth2,sth2,u_t;
  double g_tt,g_tp,g_pp;
  double gcov[NDIM][NDIM];
  double Ucon[NDIM],Ucon_bl[NDIM]; 
  double m1oA2,AA;
  double trans[NDIM][NDIM], tmp[NDIM];
  
  bl_coord(X, &r, &th);
  gcov_func(X, gcov);
  omega=1./(pow(r,1.5)+a);
  m1oA2=gcov[0][0] + 2.*omega*gcov[0][3] + omega*omega*gcov[3][3];
  AA=sqrt(-1./m1oA2);

  /* Keplerian four-velocity around Kerr black hole in BL coordinates*/
  Ucon[0] = AA;
  Ucon[1] = 0.;
  Ucon[2] = 0.;
  Ucon[3] = AA*omega;

  //gcov and gcon for MBL coordinates
  
  lower(Ucon,gcov,Ucov);
  return ;

  
}

//this is in BL coordinates
void get_model_ucon(double X[NDIM], double Ucon[NDIM])
{

  int j,k;
  double omega,l,r,r2,th,r32,r12,rho2,cth2,sth2,u_t;
  double g_tt,g_tp,g_pp;
  double gcov[NDIM][NDIM];
  double Ucon_bl[NDIM]; 
  double m1oA2,AA;
  double trans[NDIM][NDIM], tmp[NDIM];

  bl_coord(X, &r, &th);
  gcov_func(X, gcov);
  omega=1./(pow(r,1.5)+a);
  m1oA2=gcov[0][0] + 2.*omega*gcov[0][3] + omega*omega*gcov[3][3];
  AA=sqrt(-1./m1oA2);

  /* Keplerian four-velocity around Kerr black hole in BL coordinates*/
  Ucon[0] = AA;
  Ucon[1] = 0.;
  Ucon[2] = 0.;
  Ucon[3] = AA*omega;
  return ;
	
}


void get_model_bcov(double X[NDIM], double Bcov[NDIM])
{

    Bcov[0] = 0. ;
    Bcov[1] = 0. ;
    Bcov[2] = 0. ;
    Bcov[3] = 0. ;
    

}

/*should be toroidal field*/
void get_model_bcon(double X[NDIM],double Bcon[NDIM])
{

    Bcon[0] = 0. ;
    Bcon[1] = 0. ;
    Bcon[2] = 0. ;
    Bcon[3] = 0. ;

    return ;
		
}


double get_model_thetae(double X[NDIM])
{

    return 0.0;

}

double get_model_b(double X[NDIM])
{
    return 0.0;
}

double get_model_ne(double X[NDIM])
{
      return 0.;
}



void set_units(char *munitstr)
{
	double Munit=1;
	double MBH;
	
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
    MBH = 10.;
#endif
	
        /** input parameters appropriate to Sgr A* **/
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


double theta_func(double x[NDIM])
{
    double r,th;
    bl_coord(x, &r, &th);
    return th;
}

/* this is from illinois version that does not need second derivative of theta, works for funky cooridnates or whatever */
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

