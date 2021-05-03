#include "decs.h"

/* first geometry used in the model and then HARM model specification routines */

/*                  model-dependent
																				
	    geometry
	    routines:
	    risco_calc
	    rhorizon_calc
	    gcov
	    get_connection           
*/

//chose field geometry in the model
#define TOROIDAL (0)
#define VERTICAL (1)
#define RADIAL (0)

//Kerr metric only
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

//Kerr metric only
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
    dxdX[2][2] = M_PI  + hslope*2*M_PI*cos(2*M_PI*X[2]); 
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

void init_model(char *args[])
{
	set_units(args[4]);
	a = 0.0;
	Risco = risco_calc(1);
	Rh = rhorizon_calc(1);
	fprintf(stdout,"Risco=%g \n",Risco);
	Rout=100.;
	hslope=0.0;
	th_beg=0.0;
}

/* 

	these supply basic model data

*/

void get_model_ucon(double X[NDIM], double Ucon[NDIM])
{
  int j,k;
  double omega,l,r,th,rho2,cth2,sth2;
  double g_tt,g_tp,g_pp;
  double gcov[NDIM][NDIM];
  double Ucon_bl[NDIM];
  double m1oA2,AA;
  double trans[NDIM][NDIM], tmp[NDIM];
  double gcon[NDIM][NDIM];
  double Ucov[NDIM];
    
  bl_coord(X, &r, &th);

  //Keplerian velocity down to ISCO
  if(r>Risco){
    sth2=sin(th)*sin(th);
    cth2=cos(th)*cos(th);
    rho2=r*r+a*a*cth2;

    double cth=cos(th);
    double sth=sin(th);
    double s2=sth*sth;
    double a2 = a * a;
    double r2 = r * r;
    double DD = 1. - 2. / r + a2 / r2;
    double mu = 1. + a2 * cth * cth / r2;
    g_tt = -(1. - 2. / (r * mu));
    g_tp = -2. * a * s2 / (r * mu);
    g_pp = r2 * sth * sth * (1. + a2 / r2 + 2. * a2 * s2 / (r2 * r * mu));  
    omega=1./(pow(r,1.5)+a);
    m1oA2=g_tt + 2.*omega*g_tp + omega*omega*g_pp;
    AA=sqrt(-1./m1oA2);
    
    Ucon_bl[0] = AA;
    Ucon_bl[1] = 0.;
    Ucon_bl[2] = 0.;
    Ucon_bl[3] = AA*omega;
  
    /* transform Ucon from BL -> KS */
    double dtKS_drBL   = 2. * r / (r*r - 2.*r + a*a);
    double dphiKS_drBL = a / (r*r - 2.*r + a*a);
    Ucon[0] = Ucon_bl[0] +  Ucon_bl[1] * dtKS_drBL;
    Ucon[1] = Ucon_bl[1];
    Ucon[2] = Ucon_bl[2];
    Ucon[3] = Ucon_bl[3] +  Ucon_bl[1] * dphiKS_drBL;
    
    /* transform KS -> KS' coords, Ucon has to be in the same coord as K^mu */
    Ucon[1] /= r;
    Ucon[2] /= M_PI;
    
    if( isnan(Ucon[0]) ){
      fprintf(stdout,"here r=%g AA=%g omega=%g %g %g %g  %g\n",r,AA,omega,
	      g_tt,g_tp,g_pp,m1oA2);
      exit(1);
    }

  }else{ 

    gcov_func(X, gcov);
    gcon_func(gcov, gcon) ;

    Ucov[0] = -1./sqrt(-gcon[0][0]) ;
    Ucov[1] = 0. ;
    Ucov[2] = 0. ;
    Ucov[3] = 0. ;
    lower(Ucov, gcon, Ucon);

    if( isnan(Ucon[0]) ){
      fprintf(stdout,"here r=%g th=%g gcov_tt=%g\n",r,th,gcov[0][0]);
      exit(1);
    }
    
  }

  return;

}


void get_model_ucov(double X[NDIM], double Ucov[NDIM])
{  
  double Ucon[NDIM],gcov[NDIM][NDIM];
  get_model_ucon(X,Ucon);
  gcov_func(X, gcov) ;
  lower(Ucon, gcov, Ucov);
  return;
}

void get_model_bcon(double X[NDIM], double Bcon[NDIM])
{

  double Ucov[NDIM],Ucon[NDIM],gcov[NDIM][NDIM];
  double B_1,B_2,B_3;
  
  get_model_ucov(X,Ucov);
  get_model_ucon(X,Ucon);
  gcov_func(X, gcov);

  //TOROIDAL FIELD
  if(TOROIDAL){
    B_1=0.0;
    B_2=0.0;
    B_3=1.0;
  }


  if(VERTICAL){

    //zone centered gdet
    double g = gdet_func(gcov);
    //compute Aphi and via constrainted transport (CT) B field, divergence free B
    double r, th;
    bl_coord(X, &r, &th);
    double dx1=0.025;
    double dx2=0.025;
    
    // vertical
    double F11 = pow(exp(log(r)-dx1) * sin(th-dx2*M_PI),1);
    double F12 = pow(exp(log(r)-dx1) * sin(th+dx2*M_PI),1);
    double F21 = pow(exp(log(r)+dx1) * sin(th-dx2*M_PI),1);
    double F22 = pow(exp(log(r)+dx1) * sin(th+dx2*M_PI),1);

    /* flux-ct */
    B_1 = -( F11
	     - F12
	     + F21
	     - F22
	     )/(2.*dx2*g) ;

    B_2 = (  F11
	     + F12
	     - F21
	     - F22
	     )/(2.*dx1*g) ;

    B_3 = 0.0 ;
  }

  if(RADIAL){

    //zone centered gdet
    double g = gdet_func(gcov);
    //compute Aphi and via constrainted transport (CT) B field, divergence free B
    double r, th;
    bl_coord(X, &r, &th);
    double dx1=0.025;
    double dx2=0.025;

    //radial, monopolar
    double F11 = 1-cos(th-dx2*M_PI);
    double F12 = 1-cos(th+dx2*M_PI);
    double F21 = 1-cos(th-dx2*M_PI);
    double F22 = 1-cos(th+dx2*M_PI);
    B_1 = -( F11
	     - F12
	     + F21
	     - F22
	     )/(2.*dx2*g) ;

    B_2 = (  F11
	     + F12
	     - F21
	     - F22
	     )/(2.*dx1*g) ;
  
    B_3 = 0.0;
  }
  
  Bcon[0] = B_1 * Ucov[1] + B_2 * Ucov[2] + B_3 * Ucov[3] ;
  Bcon[1] = ( B_1 + Bcon[0] * Ucon[1]) / Ucon[0];
  Bcon[2] = ( B_2 + Bcon[0] * Ucon[2]) / Ucon[0];
  Bcon[3] = ( B_3 + Bcon[0] * Ucon[3]) / Ucon[0];

}

void get_model_bcov(double X[NDIM], double Bcov[NDIM])
{

  double Bcon[NDIM],gcov[NDIM][NDIM];
  get_model_bcon(X,Bcon);
  gcov_func(X, gcov);
  lower(Bcon, gcov, Bcov);

}

double get_model_thetae(double X[NDIM])
{
  
  double r, th;
  bl_coord(X, &r, &th);
  double thetae0=200.;
  return thetae0*pow(r,-0.84);

}

double get_model_b(double X[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);
  return 100.*pow(r,-1);
}

double get_model_ne(double X[NDIM])
{

  
  double r, th;
  double ne_bg,mu2;
  
  bl_coord(X, &r, &th);
  mu2=cos(th)*cos(th);
  //paramters
  double sigma2=0.3*0.3;
  ne_bg=7e3;

  if(r>Risco){
    return ne_bg*pow(r,-1.5)*exp(-mu2/2./sigma2) ;
  }else{
    return 1e-3;
  }
  

}

/* return boyer-lindquist coordinate of point */
void bl_coord(double *X, double *r, double *th)
{

	*r = exp(X[1]);
	*th = M_PI*X[2]  ;

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


double theta_func(double X[NDIM])
{
    double r,th;
    bl_coord(X, &r, &th);
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

