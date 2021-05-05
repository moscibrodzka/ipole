/**********************************************************/
/*** all you need to make a polarized radiative transfer***/
/******** used in ipole to evolve complex tensor N ********/
/******* along with standard evolution for I scalar *******/
/**********************************************************/
/**** written by Monika Moscibrodzka on 09 July 2014 ******/
/************  last update: 9 May 2017   ******************/
/*************** co-author: C.F. Gammie *******************/
/**********************************************************/

#include "decs.h"
#include "defs.h"
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>

/* the following definitions are used only locally */
#define S2     (1.41421356237310)	//sqrt(2)
#define S3     (1.73205080756888)	//sqrt(3)


/* declarations of local functions */
/* thermal plasma emissivity, absorptivity and Faraday conversion and rotation */
double g(double Xe);
double h(double Xe);
double Je(double Xe);
double jffunc(double Xe);
double I_I(double x);
double I_Q(double x);
double I_V(double x);
double besselk_asym(int n, double x);


/*************************SUPPORTING FUNCTIONS******************************/

/*invariant plasma emissivities/abs/rot in tetrad frame */
void jar_calc(double X[NDIM], double Kcon[NDIM],
	      double *jI, double *jQ, double *jU, double *jV,
	      double *aI, double *aQ, double *aU, double *aV,
	      double *rQ, double *rU, double *rV)
{
    double nu, Thetae, Ne, B, theta, nusq;
    double x, Xe, omega0, nuc;
    double Bnuinv;
    double Ucov[NDIM];
    double Thetaer, wp2;

    Ne = get_model_ne(X);
    get_model_ucov(X, Ucov);
    theta = get_bk_angle(X, Kcon, Ucov);	/* angle between k & b  */
    B = get_model_b(X);	/* field in G                */
    double sigma_m=(B*B/B_unit/B_unit)/(Ne*(MP+ME))*RHO_unit;
    nu = get_fluid_nu(Kcon, Ucov);	/* freqcgs1;  freq in Hz */
    
    if (Ne <= 0.) {  

	 *jI = 0.0;
	 *jQ = 0.0;
	 *jU = 0.0;
	 *jV = 0.0;
	 
	 *aI = 0.0;
	 *aQ = 0.0;
	 *aU = 0.0;
	 *aV = 0.0;
	 
	 *rQ = 0.0;
	 *rU = 0.0;
	 *rV = 0.0;

	 return;
     } else if (theta <= 0. || theta >= M_PI) {	/* no emission/absorption along field  */

	*jI = 0.0;
	*jQ = 0.0;
	*jU = 0.0;
	*jV = 0.0;

	*aI = 0.0;
	*aQ = 0.0;
	*aU = 0.0;
	*aV = 0.0;

	*rQ = 0.0;
	*rU = 0.0;

	nu = get_fluid_nu(Kcon, Ucov);	/* freqcgs1;  freq in Hz */
	wp2 = 4. * M_PI * Ne * EE * EE / ME;
	B = get_model_b(X);		/* field in G                */
	omega0 = EE * B / ME / CL;
	Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
	Thetaer = 1. / Thetae;
					/* Faraday rotativities for thermal plasma */
	Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));

	// other possibilities
        //original Dexter 2016 
        //*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	//  (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
	//with my correction
	//*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	//    (gsl_sf_bessel_Kn(0,1./Thetae) -Je(Xe)) / gsl_sf_bessel_Kn(2,1./Thetae) * cos(theta);
	//Hung and Scherbakov 2011 fit
	//*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	//   gsl_sf_bessel_Kn(0, Thetaer) / (gsl_sf_bessel_Kn(2, Thetaer)+SMALL) * g(Xe) * cos(theta);
	
	//this is used for EHT-library only
	if (Thetae > 3.0) {
	    // High temperature: use approximations to bessel
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
		(besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
	} else if (0.2 < Thetae && Thetae <= 3.0) {
	    // Mid temperature: use real bessel functions 
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
		(gsl_sf_bessel_Kn(0, Thetaer) - Je(Xe)) / gsl_sf_bessel_Kn(2, Thetaer) * cos(theta);
	} else if (Thetae <= 0.2) {
	    // Use the constant low-temperature limit
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) * cos(theta);
	}
		
	*rV *= nu;
	 
	return;
	
    } else {

      nu = get_fluid_nu(Kcon, Ucov);	/* freqcgs1;  freq in Hz */
      nusq = nu * nu;
      B = get_model_b(X);	/* field in G                */
      Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
      Thetaer = 1. / Thetae;

      omega0 = EE * B / ME / CL;
      wp2 = 4. * M_PI * Ne * EE * EE / ME;

      /* Faraday rotativities for thermal plasma */
      Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));

      /* Here I use approximate bessel functions to match rhoqv with grtrans */
      *rQ = 2. * M_PI * nu / 2. / CL * wp2 * omega0 * omega0 / pow(2 * M_PI * nu, 4) *
	jffunc(Xe) * (besselk_asym(1, Thetaer) / besselk_asym(2, Thetaer) +
		      6. * Thetae) * sin(theta) * sin(theta);
      
      *rU = 0.0;

      //original Dexter 2016 
      //*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      //	(besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
      //my correction
      //*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      //	  (gsl_sf_bessel_Kn(0,1./Thetae) - Je(Xe)) / gsl_sf_bessel_Kn(2,1./Thetae) * cos(theta);
      //Hung and Scherbakov 2011 fit
      //*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      //   gsl_sf_bessel_Kn(0, Thetaer) / (gsl_sf_bessel_Kn(2, Thetaer)+SMALL) * g(Xe) * cos(theta);  
      
      //this is used for EHT-library only
      if (Thetae > 3.0) {
	// High temperature: use approximations to bessel
	*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	  (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
      } else if (0.2 < Thetae && Thetae <= 3.0) {
	// Mid temperature: use real bessel functions (TODO fit?)
	*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	  (gsl_sf_bessel_Kn(0, Thetaer) - Je(Xe)) / gsl_sf_bessel_Kn(2, Thetaer) * cos(theta);
      } else if (Thetae <= 0.2) {
	// Use the constant low-temperature limit
	*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) * cos(theta);
      }
      

      /* invariant rotativities */
      *rQ *= nu;
      *rV *= nu;
      
      /* synchrotron emissivity */ 
      nuc = 3.0 * EE * B * sin(theta) / 4.0 / M_PI / ME / CL * Thetae * Thetae + 1.0;
      x = nu / nuc;
      
      *jI = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_I(x); // [g/s^2/cm = ergs/s/cm^3]
      *jQ = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_Q(x);
      *jU = 0.0;	                                        	 // convention; depends on tetrad
      *jV = 2. * Ne * EE * EE * nu / tan(theta) / 3. / S3 / CL / Thetae / Thetae / Thetae * I_V(x);

      /* invariant emissivity */
      *jI /= nusq;
      *jQ /= nusq;
      *jU /= nusq;
      *jV /= nusq;
      
      /* invariant synchrotron absorptivity */
      Bnuinv = Bnu_inv(nu, Thetae) + 1e-100;   /* Planck function */
      *aI = *jI / Bnuinv;
      *aQ = *jQ / Bnuinv;
      *aU = *jU / Bnuinv;
      *aV = *jV / Bnuinv;
            
      return;
      
    }



}

/*invariant plasma emissivities/abs/rot in tetrad frame */
/* for mixed thermal+power-law eDF*/
void jar_calc_mixed_pl(double X[NDIM], double Kcon[NDIM],
	      double *jI, double *jQ, double *jU, double *jV,
	      double *aI, double *aQ, double *aU, double *aV,
	      double *rQ, double *rU, double *rV)
{
    double nu, Thetae, Ne, B, theta, nusq;
    double x, Xe, omega0, nuc;
    double Bnuinv;
    double Ucov[NDIM];
    double Thetaer, wp2;

    Ne = get_model_ne(X);
    get_model_ucov(X, Ucov);
    theta = get_bk_angle(X, Kcon, Ucov);	/* angle between k & b  */
    
    if (Ne <= 0.) {  // avoid 1./0. issues

	 *jI = 0.0;
	 *jQ = 0.0;
	 *jU = 0.0;
	 *jV = 0.0;
	 
	 *aI = 0.0;
	 *aQ = 0.0;
	 *aU = 0.0;
	 *aV = 0.0;
	 
	 *rQ = 0.0;
	 *rU = 0.0;
	 *rV = 0.0;

	 return;
	 
     } else if (theta <= 0. || theta >= M_PI) {	/* no emission/absorption along field  */

	*jI = 0.0;
	*jQ = 0.0;
	*jU = 0.0;
	*jV = 0.0;

	*aI = 0.0;
	*aQ = 0.0;
	*aU = 0.0;
	*aV = 0.0;

	*rQ = 0.0;
	*rU = 0.0;

	nu = get_fluid_nu(Kcon, Ucov);	/* freqcgs1;  freq in Hz */
	wp2 = 4. * M_PI * Ne * EE * EE / ME;
	B = get_model_b(X);		/* field in G                */
	omega0 = EE * B / ME / CL;
	Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
	Thetaer = 1. / Thetae;
					/* Faraday rotativities for thermal plasma */
	Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));


	//original Dexter 2016
	//*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	//        (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
	//my correction
	//*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	//          (gsl_sf_bessel_Kn(0,1./Thetae) - Je(Xe)) / gsl_sf_bessel_Kn(2,1./Thetae) * cos(theta);
	//Hung and Scherbakov 2011 fit
	//*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	//   gsl_sf_bessel_Kn(0, Thetaer) / (gsl_sf_bessel_Kn(2, Thetaer)+SMALL) * g(Xe) * cos(theta);

	//this is used for EHT-library only
	if (Thetae > 3.0) {
	    // High temperature: use approximations to bessel
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
		(besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
	} else if (0.2 < Thetae && Thetae <= 3.0) {
	    // Mid temperature: use real bessel functions 
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
		(gsl_sf_bessel_Kn(0, Thetaer) - Je(Xe)) / gsl_sf_bessel_Kn(2, Thetaer) * cos(theta);
	} else if (Thetae <= 0.2) {
	    // Use the constant low-temperature limit
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) * cos(theta);
	}


       /* invariant rotativities */
	*rV *= nu;

	return;

     } else {

       nu = get_fluid_nu(Kcon, Ucov);	/* freqcgs1;  freq in Hz */
       nusq = nu * nu;
       B = get_model_b(X);	/* field in G                */
       Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
       Thetaer = 1. / Thetae;

       omega0 = EE * B / ME / CL;
       wp2 = 4. * M_PI * Ne * EE * EE / ME;

       /* Faraday rotativities for thermal plasma */
       Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));


       
      /* Here I use approximate bessel functions to match rhoqv with grtrans */
      *rQ = 2. * M_PI * nu / 2. / CL * wp2 * omega0 * omega0 / pow(2 * M_PI * nu, 4) *
	jffunc(Xe) * (besselk_asym(1, Thetaer) / besselk_asym(2, Thetaer) +
		      6. * Thetae) * sin(theta) * sin(theta);
      
      *rU = 0.0;

      //original Dexter 2016
      //*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      //        (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
      //my correction
      //*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      //          (gsl_sf_bessel_Kn(0,1./Thetae) - Je(Xe)) / gsl_sf_bessel_Kn(2,1./Thetae) * cos(theta);
      //Hung and Scherbakov 2011 fit
      //*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      //   gsl_sf_bessel_Kn(0, Thetaer) / (gsl_sf_bessel_Kn(2, Thetaer)+SMALL) * g(Xe) * cos(theta);
      
      //this is used for EHT-library only
      if (Thetae > 3.0) {
	  // High temperature: use approximations to bessel
	  *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	      (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
      } else if (0.2 < Thetae && Thetae <= 3.0) {
	  // Mid temperature: use real bessel functions 
	  *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	      (gsl_sf_bessel_Kn(0, Thetaer) - Je(Xe)) / gsl_sf_bessel_Kn(2, Thetaer) * cos(theta);
      } else if (Thetae <= 0.2) {
	  // Use the constant low-temperature limit
	  *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) * cos(theta);
      }
	
      /* invariant rotativities */
      *rQ *= nu;
      *rV *= nu;

      /*synchrotron emissivity */
      nuc = 3.0 * EE * B * sin(theta) / 4.0 / M_PI / ME / CL * Thetae * Thetae + 1.0;
      x = nu / nuc;

      *jI = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_I(x); // [g/s^2/cm = ergs/s/cm^3]
      *jQ = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_Q(x); 
      *jU = 0.0;	                                        	 // convention; depends on tetrad
      *jV = 2. * Ne * EE * EE * nu / tan(theta) / 3. / S3 / CL / Thetae / Thetae / Thetae * I_V(x);
      
      /* invariant emissivity */
      *jI /= nusq;
      *jQ /= nusq;
      *jU /= nusq;
      *jV /= nusq;
      
      /* invariant synchrotron absorptivity */
      Bnuinv = Bnu_inv(nu, Thetae)+1e-100;   /* Planck function */
      *aI = *jI / Bnuinv;
      *aQ = *jQ / Bnuinv;
      *aU = *jU / Bnuinv;
      *aV = *jV / Bnuinv;

      /////// NONTHERMAL TAIL /////////
      /* add some nonthermal component to thermal to form hybrid, gmin should be computed for given p*/
      double p=3.5;
      double eta=0.1;
      double gmin=9.*Thetae; //this is a paramter that has to be evaluated
      double gmax=1e6;
      
      //add nonthermal emissivity, abs, adjust npl
      nuc = EE*B/(2.*M_PI*ME*CL);
      double a_thetae=(6. + 15.*Thetae)/(4. + 5.*Thetae);
      double Ne_pl=(p-2.)/(p-1.)/gmin*eta*a_thetae*Thetae*Ne;
      double emcon=Ne_pl*EE*EE*nuc/CL;
      double j_pl,j_pl_I,j_pl_Q,j_pl_U,j_pl_V;

      j_pl = emcon*pow(3,p/2.)*(p-1)*sin(theta)/2./(p+1)/(pow(gmin,1-p)-pow(gmax,1-p))*
	gsl_sf_gamma((3*p-1)/12)*gsl_sf_gamma((3*p+19)/12)*
	pow(nu/nuc/sin(theta),-(p-1)/2);
       
      j_pl_I = j_pl;
      j_pl_Q = j_pl*(p+1.)/(p+7./3.);
      j_pl_U = 0.0;
      j_pl_V = 171./250.*sqrt(p)/tan(theta)*pow(nu/3./nuc/sin(theta),-0.5)*j_pl;

      /*add invariant emiss to thermal DF*/
      *jI += j_pl_I/nusq;
      *jQ += j_pl_Q/nusq;
      *jU += j_pl_U/nusq;
      *jV += j_pl_V/nusq;
      
       /* absorption */
      double a_pl,a_pl_I,a_pl_Q,a_pl_U,a_pl_V;
      //there was a bug but now it is fixed
      emcon=Ne_pl*EE*EE/nu/ME/CL;
       
      a_pl=emcon*pow(3,(p+1)/2.)*(p-1)/4./(pow(gmin,1-p)-pow(gmax,1-p))*
	gsl_sf_gamma((3*p+2)/12)*gsl_sf_gamma((3*p+22)/12)*
	pow(nu/nuc/sin(theta),-(p+2)/2);
      a_pl_I=a_pl;
      a_pl_Q=a_pl*0.75*pow(p-1,43/500);
      a_pl_U=0.0;
      a_pl_V=7/4.*pow(71/100*p+22/625,197/500)*
	pow(pow(sin(theta),-48/25)-1,64/125)*
	pow(nu/nuc/sin(theta),-0.5)*a_pl;

      /*add invariant abs to thermal DF*/
      *aI += a_pl_I*nu;
      *aQ += a_pl_Q*nu;
      *aU += a_pl_U*nu;
      *aV += a_pl_V*nu;
      
      /*add rotatativity for power-law*/
      double nub=EE*B/2./M_PI/ME/CL;
      double Cv=2.*(p+2.)/(p+1.);
      double rho1=Ne_pl*EE*EE/ME/CL/nub/sin(theta)*(p-1.)/(pow(gmin,1.-p)-pow(gmax,1.-p));
      double numin=gmin*gmin*nub*sin(theta);//ok
      
      /* add invariant rotativity to thermal component, here sign of rQ is fixed w/ respect to Dexter 2016, Johnson 977*/
      *rQ = *rQ + (rho1 * pow(nub*sin(theta)/nu,3) * pow(gmin,2-p) * (1.-pow(numin/nu,p/2.-1.)) *pow(p/2.-1.,-1)) *nu ;
      *rV = *rV + (Cv *rho1*pow(nub*sin(theta)/nu,2)*pow(gmin,-p-1)*log(gmin)/tan(theta)) * nu;
    
    }
      
}



/* invariant plasma emissivities/abs/rot in tetrad frame */
/* for thermal and kappa distribution */
/* adding all emissivities, abs and rotativities*/
void jar_calc_mixed_kappa(double X[NDIM], double Kcon[NDIM],
	      double *jI, double *jQ, double *jU, double *jV,
	      double *aI, double *aQ, double *aU, double *aV,
	      double *rQ, double *rU, double *rV)
{
    double nu, Thetae, Ne, B, theta, nusq;
    double x, Xe, omega0, nuc;
    double Bnuinv;
    double Ucov[NDIM];
    double Thetaer, wp2;

    Ne = get_model_ne(X);
    get_model_ucov(X, Ucov);
    theta = get_bk_angle(X, Kcon, Ucov);	/* angle between k & b  */

    /* paramters of KAPPA function*/
    //double eta=0.0; //by default it is zero
    //or bernoulli paramters
    //double rho1=Ne/RHO_unit*(MP+ME); //dimentionless rho
    //    double uu=get_model_uu(X);
    //    double Be=-(1+uu/rho1*gam)*Ucov[0];
    //    if(Be > 1.02) eta=0.1; //along the jet
    //eta=0.01;//everywhere
    /* END PARAMETERS OF KAPPA function*/

    B = get_model_b(X);		/* field in G                */
    double sigma_m=(B*B/B_unit/B_unit)/(Ne*(MP+ME))*RHO_unit;
    nu = get_fluid_nu(Kcon, Ucov);	/* freqcgs1;  freq in Hz */
    
    if (Ne <= 0.) {  // avoid 1./0. issues

	 *jI = 0.0;
	 *jQ = 0.0;
	 *jU = 0.0;
	 *jV = 0.0;
	 
	 *aI = 0.0;
	 *aQ = 0.0;
	 *aU = 0.0;
	 *aV = 0.0;
	 
	 *rQ = 0.0;
	 *rU = 0.0;
	 *rV = 0.0;

	 return;
	 
    } else if (theta <= 0. || theta >= M_PI) {	/* no emission/absorption along field  */

	*jI = 0.0;
	*jQ = 0.0;
	*jU = 0.0;
	*jV = 0.0;

	*aI = 0.0;
	*aQ = 0.0;
	*aU = 0.0;
	*aV = 0.0;

	*rQ = 0.0;
	*rU = 0.0;

	nu = get_fluid_nu(Kcon, Ucov);	/* freqcgs1;  freq in Hz */
	wp2 = 4. * M_PI * Ne * EE * EE / ME;
	B = get_model_b(X);		/* field in G                */
	omega0 = EE * B / ME / CL;
	Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
	Thetaer = 1. / Thetae;
					/* Faraday rotativities for thermal plasma */
	Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));

	*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	    (gsl_sf_bessel_Kn(0,1./Thetae) -Je(Xe)) / gsl_sf_bessel_Kn(2,1./Thetae) * cos(theta);

	//used for EHT-library only
	/*
	if (Thetae > 3.0) {
	    // High temperature: use approximations to bessel
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
		(besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
	} else if (0.2 < Thetae && Thetae <= 3.0) {
	    // Mid temperature: use real bessel functions (TODO fit?)
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
		(gsl_sf_bessel_Kn(0, Thetaer) - Je(Xe)) / gsl_sf_bessel_Kn(2, Thetaer) * cos(theta);
	} else if (Thetae <= 0.2) {
	    // Use the constant low-temperature limit
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) * cos(theta);
	}
	*/
	
	/* invariant rotativities */
	*rV *= nu;
	return;
	

    } else {

      nu = get_fluid_nu(Kcon, Ucov);	/* freqcgs1;  freq in Hz */
      nusq = nu * nu;
      B = get_model_b(X);	/* field in G                */
      Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
      Thetaer = 1. / Thetae;

      // but if sigma_m>1, zero thermal coeffs and add kappa only
      double sigma_m=(B*B/B_unit/B_unit)/(Ne/RHO_unit*(MP+ME)); 
      
      // outside of the funnel -> thermal plasma emission, absorption and rotation
      if(sigma_m <= 0.85){

	omega0 = EE * B / ME / CL;
	wp2 = 4. * M_PI * Ne * EE * EE / ME;
	
	/* Faraday rotativities for thermal plasma */
	Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));
	
	/* Here I use approximate bessel functions to match rhoqv with grtrans */
	*rQ = 2. * M_PI * nu / 2. / CL * wp2 * omega0 * omega0 / pow(2 * M_PI * nu, 4) *
	  jffunc(Xe) * (besselk_asym(1, Thetaer) / besselk_asym(2, Thetaer) +
			6. * Thetae) * sin(theta) * sin(theta);
	*rU = 0.0;

	//my stuff
	//      *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	//       (gsl_sf_bessel_Kn(0,1./Thetae) - Je(Xe)) / gsl_sf_bessel_Kn(2,1./Thetae) * cos(theta);
	
	// Switch between three different fits for rho_V
	//EHT-library
	if (Thetae > 3.0) {
	  // High temperature: use approximations to bessel
	  *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	    (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
	} else if (0.2 < Thetae && Thetae <= 3.0) {
	  // Mid temperature: use real bessel functions (TODO fit?)
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
		(gsl_sf_bessel_Kn(0, Thetaer) - Je(Xe)) / gsl_sf_bessel_Kn(2, Thetaer) * cos(theta);
	} else if (Thetae <= 0.2) {
	    // Use the constant low-temperature limit
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) * cos(theta);
	}

	/* invariant rotativities */
	*rQ *= nu;
	*rU *= nu;
	*rV *= nu;

	/*synchrotron emissivity */
	nuc = 3.0 * EE * B * sin(theta) / 4.0 / M_PI / ME / CL * Thetae * Thetae + 1.0;
	x = nu / nuc;
	*jI = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_I(x);	// [g/s^2/cm = ergs/s/cm^3]
	*jQ = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_Q(x); //here mistake in grtrans
	*jU = 0.0;		// convention; depends on tetrad
	*jV = 2. * Ne * EE * EE * nu / tan(theta) / 3. / S3 / CL / Thetae / Thetae / Thetae * I_V(x);
         
	/* invariant emissivity */
	*jI /= nusq;
	*jQ /= nusq;
	*jU /= nusq;
	*jV /= nusq;
      
	/* invariant synchrotron absorptivity */
	Bnuinv = Bnu_inv(nu, Thetae) + 1e-100;   /* Planck function */
	*aI = *jI / Bnuinv;
	*aQ = *jQ / Bnuinv;
	*aU = *jU / Bnuinv;
	*aV = *jV / Bnuinv;
	
	//regions in sigma_m > 1, within the jet, only emission,abs
      } else {

	*jI = 0.0;
	*jQ = 0.0;
	*jU = 0.0;
	*jV = 0.0;

	*aI = 0.0;
	*aQ = 0.0;
	*aU = 0.0;
	*aV = 0.0;
	
	*rQ = 0.0;
	*rU = 0.0;
	*rV = 0.0;

	//comment just this out because below we are adding
	//       	return;

	/********************* EMISSION, ABS only ******************/
	/*    /\* add eta nonthermal particles in kappa distribution *\/ */
	
	double kappa=3.5;  
	//width of the distribution function
	//either from thermal
	//double w=Thetae;
	
	//or a fraction of magnetic energy
	double eta=1.0; 
	double w=eta*(kappa-3.)/6./kappa*MP/ME*sigma_m;
	  
	//stokes I
	//nuc changes meaning
	double nu_c=EE*B/(2.*M_PI*ME*CL) ;
	double nu_w = nu_c*pow(w*kappa,2.)*sin(theta);
	double X_k = nu/nu_w;
	double prefactor = (Ne* pow(EE, 2.) * nu_c
			    * sin(theta))/CL;
	double Nlow = 4. * M_PI * tgamma(kappa-4./3.)
	  / (pow(3., 7./3.) * tgamma(kappa-2.));
	double Nhigh = (1./4.) * pow(3., (kappa-1.)/2.)
	  * (kappa-2.) * (kappa-1.)
	  * tgamma(kappa/4.-1./3.)
	  * tgamma(kappa/4.+4./3.);
	double x = 3. * pow(kappa, -3./2.);
	double ans = prefactor * Nlow * pow(X_k, 1./3.)
	  * pow(1.+pow(X_k, x * (3. * kappa-4.)/6.)
		* pow(Nlow/Nhigh, x), -1./x);

	//invariant !
	*jI += ans/nu/nu;
      
	//stokes Q
	Nlow = -(1./2.) * 4. * M_PI * tgamma(kappa-4./3.)
	  /(pow(3., 7./3.) * tgamma(kappa-2.));
	Nhigh = -(pow(4./5., 2)+kappa/50.) * (1./4.)
	    * pow(3., (kappa-1.)/2.) * (kappa-2.)
	    * (kappa-1.) * tgamma(kappa/4.-1./3.)
	    * tgamma(kappa/4.+4./3.);
	  x = (37./10.)*pow(kappa, -8./5.);
	  ans = prefactor * Nlow * pow(X_k, 1./3.)
	    * pow(1. + pow(X_k, x * (3. * kappa-4.)/6.)
		  * pow(Nlow/Nhigh, x), -1./x);
	  
	  //   fprintf(stdout,"theta=%g jq= %g %g \n",Thetae, *jQ, eta*ans/nu/nu);
	  //here I change sign
	  //invariant 
	  *jQ += -ans/nu/nu;

	  //stokes U
	  *jU += 0.0;
	  
	  //stokes V
	  Nlow = -pow(3./4., 2.)
	    * pow(pow(sin(theta), -12./5.)-1., 12./25.)
	    * (pow(kappa, -66./125.) / w)
	    * pow(X_k, -7./20.) * 4. * M_PI
	    * tgamma(kappa-4./3.)
	    / (pow(3., 7./3.)*tgamma(kappa-2.));
      
	  Nhigh = -pow(7./8., 2.)
	    * pow(pow(sin(theta), -5./2.)-1., 11./25.)
	    * (pow(kappa, -11./25.) / w)
	    * pow(X_k, -1./2.) * (1./4.) * pow(3., (kappa-1.)/2.)
	    * (kappa-2.) * (kappa-1.)
	    * tgamma(kappa/4.-1./3.)
	    * tgamma(kappa/4.+4./3.);
      
	  x = 3.*pow(kappa, -3./2.);

	  ans = prefactor * Nlow * pow(X_k, 1./3.)
	    * pow(1.+pow(X_k, x * (3.*kappa-4.)/6.)
		  * pow(Nlow/Nhigh, x), -1./x);
      
	  // fprintf(stdout,"jv= %g %g \n", *jV, eta*ans/nu/nu);
	  //invariant
	  *jV += ans/nu/nu;
	  //end emissivities, kappa emissivities added


	  
	  //stokes I absorption
	  prefactor = Ne * EE  / (B * sin(theta));
	  //this cannot be...
	  double aaa = kappa - 1./3.;
	  double b = kappa + 1.;
	  double c = kappa + 2./3.;
	  double z = -kappa*w;
      
	  /* GSL 2F1 only works for |z| < 1; had to apply a hypergeometric function
	     identity because in our case z = -kappa*w, so |z| > 1 */
	  
	  double hyp2f1 = pow(1.-z, -aaa) * tgamma(c) * tgamma(b-aaa)
	    / (tgamma(b)*tgamma(c-aaa))
	    * gsl_sf_hyperg_2F1(aaa, c-b, aaa-b+1., 1./(1.-z))
	    + pow(1.-z, -b) * tgamma(c) * tgamma(aaa-b)
	    / (tgamma(aaa) * tgamma(c-b))
	    * gsl_sf_hyperg_2F1(b, c-aaa, b-aaa+1., 1./(1.-z));
	  
	  Nlow = pow(3., 1./6.) * (10./41.) * pow(2. * M_PI, 2.)
	    / pow(w * kappa, 16./3.-kappa)
	    * (kappa-2.) * (kappa-1.) * kappa
	    / (3.*kappa-1.) * tgamma(5./3.) * hyp2f1;
	  
	  Nhigh = 2. * pow(M_PI, 5./2.)/3. * (kappa-2.)
	    * (kappa-1.) * kappa
	    / pow(w * kappa, 5.)
	    * (2 * tgamma(2. + kappa/2.)
	       / (2.+kappa)-1.)
	    * (pow(3./kappa, 19./4.) + 3./5.);
    
	  x = pow(-7./4. + 8. * kappa/5., -43./50.);
      
	  ans = prefactor * Nlow * pow(X_k, -5./3.)
	    * pow(1. + pow(X_k, x * (3. * kappa-1.)/6.)
		  * pow(Nlow/Nhigh, x), -1./x);

	  
	  *aI += ans*nu;
	  
	  //stokes Q
	  /*GSL 2F1 only works for |z| < 1; had to apply a hypergeometric function
	    identity because in our case z = -kappa*w, so |z| > 1 */
	  hyp2f1 = pow(1.-z, -aaa) * tgamma(c) * tgamma(b-aaa)
	    / (tgamma(b) * tgamma(c-aaa))
	    * gsl_sf_hyperg_2F1(aaa, c-b, aaa-b+1., 1./(1.-z))+pow(1.-z, -b)
	    * tgamma(c) * tgamma(aaa-b) / (tgamma(aaa)*tgamma(c-b))
	    * gsl_sf_hyperg_2F1(b, c - aaa, b - aaa + 1., 1./(1. - z));
	  Nlow = -(25./48.) * pow(3., 1./6.) * (10./41.)
	    * pow(2. * M_PI, 2.)
	    / pow(w*kappa, 16./3.-kappa)
	    * (kappa-2.) * (kappa-1.)
	    * kappa / (3. * kappa-1.)
	    * tgamma(5./3.) * hyp2f1;
	  Nhigh = -(pow(21., 2.) * pow(kappa, -144./25.) + 11./20.)
	    * 2. * pow(M_PI, 5./2.)/3. * (kappa-2.)
	    * (kappa-1.) * kappa
	    / pow(w * kappa, 5.)
	    * (2 * tgamma(2. + kappa/2.)
	       / (2. + kappa)-1.);
      
	  x = (7./5.) * pow(kappa, -23./20.);
	  
	  ans = prefactor * Nlow * pow(X_k, -5./3.)
	    * pow(1. + pow(X_k, x * (3. * kappa-1.) / 6.)
		  * pow(Nlow/Nhigh, x), -1./x);

	  //change of sign to match ipole conventions
	  //invariant
	  *aQ += -ans*nu;


	  
	  //stokes U
	  *aU += 0.0;
	  
	  //stokes V

	  /*GSL 2F1 only works for |z| < 1; had to apply a hypergeometric function
	    identity because in our case z = -kappa*w, so |z| > 1 */
	  hyp2f1 = pow(1.-z, -aaa) * tgamma(c) * tgamma(b-aaa)
	    / (tgamma(b) * tgamma(c-aaa))
	    * gsl_sf_hyperg_2F1(aaa, c-b, aaa-b+1., 1./(1.-z))
	    + pow(1.-z, -b) * tgamma(c) * tgamma(aaa-b)
	    / (tgamma(aaa) * tgamma(c-b))
	    * gsl_sf_hyperg_2F1(b, c - aaa, b - aaa + 1., 1./(1.-z));
      
	  Nlow = -(77./(100. * w))
	    * pow(pow(sin(theta), -114./50.)
		  -1., 223./500.) * pow(X_k, -7./20.)
	    * pow(kappa, -7./10) * pow(3., 1./6.)
	    * (10./41.) * pow(2. * M_PI, 2.)
	    / pow(w*kappa, 16./3.-kappa)
	    * (kappa-2.) * (kappa-1.)
	    * kappa / (3. * kappa-1.)
	    * tgamma(5./3.) * hyp2f1;
	  Nhigh = -(143./10. * pow(w, -116./125.))
	    * pow(pow(sin(theta), -41./20.)-1., 1./2.)
	    * (13.*13. * pow(kappa, -8.) + 13./(2500.)
	       * kappa - 263./5000. + 47.
	       / (200.*kappa)) * pow(X_k, -1./2.)
	    * 2. * pow(M_PI, 5./2.) / 3.
	    * (kappa-2.) * (kappa-1.) * kappa
	    / pow(w*kappa, 5.)
	    * (2 * tgamma(2. + kappa/2.)
	       / (2. + kappa) - 1.);
      
	  x = (61./50.)*pow(kappa, -142./125.)+7./1000.;
	  
	  ans = prefactor * Nlow * pow(X_k, -5./3.)
	    * pow(1.+pow(X_k, x * (3. * kappa-1.) / 6.)
		  * pow(Nlow/Nhigh, x), -1./x);
      
	  /*The Stokes V absorption coefficient changes sign at observer_angle
	    equals 90deg, but this formula does not.  This discrepancy is a
	    bug in this formula, and is patched by the term below.*/
	  double sign_bug_patch = cos(theta)/fabs(cos(theta));

	  /*NOTE: Sign corrected; the sign in Leung et al. (2011)
	    and Pandya et al. (2016) for Stokes V transfer coefficients
	    does not follow the convention the papers describe (IEEE/IAU);
	    the sign has been corrected here.*/
	  ans=-ans * sign_bug_patch;
	  //invariant !
	  *aV += ans*nu;
	  
      }
      
    }


}



/*emissivity functions and functions used for Faraday conversion and rotation*/
/*from J. Dexter PhD thesis (checked with Leung harmony program, and Huang & Shcherbakov 2011*/
double g(double Xe)
{
    return 1. - 0.11 * log(1 + 0.035 * Xe);
}

double h(double Xe)
{
    return 2.011 * exp(-pow(Xe, 1.035) / 4.7) -
	cos(Xe * 0.5) * exp(-pow(Xe, 1.2) / 2.73) -
	0.011 * exp(-Xe / 47.2);
}

double Je(double Xe)
{
    return 0.43793091 * log(1. + 0.00185777 * pow(Xe, 1.50316886));
}

double jffunc(double Xe)
{
    double extraterm;
    extraterm =
	(0.011 * exp(-Xe / 47.2) -
	 pow(2., -1. / 3.) / pow(3.,
				 23. / 6.) * M_PI * 1e4 * pow(Xe + 1e-16,
							      -8. / 3.)) *
	(0.5 + 0.5 * tanh((log(Xe) - log(120.)) / 0.1));
    return 2.011 * exp(-pow(Xe, 1.035) / 4.7) -
	cos(Xe * 0.5) * exp(-pow(Xe, 1.2) / 2.73) -
	0.011 * exp(-Xe / 47.2) + extraterm;
}

double I_I(double x)
{
    return 2.5651 * (1 + 1.92 * pow(x, -1. / 3.) +
		     0.9977 * pow(x, -2. / 3.)) * exp(-1.8899 * pow(x,
								    1. /
								    3.));
}

double I_Q(double x)
{
    return 2.5651 * (1 + 0.93193 * pow(x, -1. / 3.) +
		     0.499873 * pow(x, -2. / 3.)) * exp(-1.8899 * pow(x,
								      1. /
								      3.));
}

double I_V(double x)
{
    return (1.81348 / x + 3.42319 * pow(x, -2. / 3.) +
	    0.0292545 * pow(x, -0.5) + 2.03773 * pow(x,
						     -1. / 3.)) *
	exp(-1.8899 * pow(x, 1. / 3.));
}

double besselk_asym(int n, double x)
{

    if (n == 0)
	return -log(x / 2.) - 0.5772;

    if (n == 1)
	return 1. / x;

    if (n == 2)
	return 2. / x / x;

    fprintf(stderr,"this cannot happen\n");
    exit(1);
}

/* end of emissivity functions */

#undef S2
#undef S3

int radiating_region(double X[4])
{
    double ne=get_model_ne(X);
#if(NT_PROB)
    return 0;
#endif
    if(ne >0.0 && X[1]<log(1000.) && X[2] > th_beg/M_PI && X[2] < (1.-th_beg/M_PI) ) return(1);
//    if(ne >0.0 && X[1]<log(100.)  ) return(1);
    else return(0);
}


double jnu_synch(double nu, double Ne, double Thetae, double B, double theta)
{
    double K2,nuc,nus,x,f,j,sth ;

    
      if (Thetae < THETAE_MIN){
          return 0.;
      }
      if (Thetae > THETAE_MAX){
	K2=2. * Thetae * Thetae;
      }else{
	K2 = gsl_sf_bessel_Kn(2,1./Thetae);
      }
    
      //original
      //K2 = gsl_sf_bessel_Kn(2,1./Thetae);
      //K2 = 2.*Thetae*Thetae ;
      
      nuc = EE*B/(2.*M_PI*ME*CL);
      sth = sin(theta);
      nus = (2./9.)*nuc*Thetae*Thetae*sth;
      if (nu > 1.e12 * nus) return (0.);
      x = nu / nus;
      double xx,xp1;
      xp1 = pow(x, 1. / 3.);
      xx = sqrt(x) + pow(2.,11./12.) * sqrt(xp1);
      f = xx * xx;
      j = (M_SQRT2 * M_PI * EE * EE * Ne * nus / (3. * CL * K2)) * f *
	  exp(-xp1);
    
      return(j) ;

}


//from symphony code
double jnu_synchSI(double nu, double Ne, double Thetae, double B,
		   double theta)
{

  double nu_c = EE * B / (2. * M_PI * ME * CL);
  double nu_s = (2./9.)*nu_c*sin(theta)*Thetae*Thetae;
  double X = nu/nu_s;
  double prefactor = (Ne * EE*EE * nu_c)/CL;
  double term1 = sqrt(2.)*M_PI/27. * sin(theta);
  double term2 = pow(pow(X, 0.5)+pow(2., 11./12.)*pow(X, 1./6.), 2.);
  double term3 = exp(-pow(X, 1./3.));
  double ans = prefactor * term1 * term2 * term3;

  return ans;
}

//from symphony
double jnu_synchSQ(double nu, double Ne, double Thetae, double B,
		   double theta){
    
  double nu_c = EE * B / (2. * M_PI * ME * CL);
  double nu_s = (2./9.)*nu_c*sin(theta)*Thetae*Thetae;
  double X = nu/nu_s;
  double prefactor = (Ne * EE*EE * nu_c)/CL;
  double term1 = sqrt(2.)*M_PI/27. * sin(theta);
  double term2 = (7.*pow(Thetae, 24./25.)+35.)
                /(10.*pow(Thetae, 24./25.)+75.);
  double term3 = pow(pow(X, 0.5)+term2*pow(2., 11./12.)*pow(X, 1./6.), 2.);
  double ans = prefactor*term1*term3*exp(-pow(X, 1./3.));
  //  return -ans;
  return ans;

}

//from symphony
double jnu_synchSV(double nu, double Ne, double Thetae, double B,
		   double theta)
{

  double nu_c = EE * B / (2. * M_PI * ME * CL);
  double nu_s = (2./9.)*nu_c*sin(theta)*Thetae*Thetae;
  double X = nu/nu_s;
  double prefactor = (Ne * EE*EE * nu_c)/CL;
  double term1 = (37.-87.*sin(theta-28./25.))
    /(100.*(Thetae+1.));
  double term2 = pow(1.+(pow(Thetae, 3./5.)/25.+7./10.)
                     *pow(X, 9./25.), 5./3.);
  double ans = prefactor*term1*term2*exp(-pow(X, 1./3.));

  /*NOTE: Sign corrected; the sign in Leung et al. (2011)
    and Pandya et al. (2016) for Stokes V transfer coefficients
    does not follow the convention the papers describe (IEEE/IAU);
    the sign has been corrected here.*/
  return ans;
}


/* get the invariant emissivity and opacity at a given position for a given wavevector */
/* here only thermal plasma emissivities */
void get_jkinv(double X[NDIM], double Kcon[NDIM], double *jnuinv,
	       double *knuinv, double col_v[3])
{
    double nu, theta, B, Thetae, Ne, Bnuinv;
    double Ucov[NDIM], Bcov[NDIM];
    double Ucon[NDIM], Bcon[NDIM];
    double Kcov[NDIM], gcov[NDIM][NDIM];

    /* get fluid parameters */
    Ne = get_model_ne(X);	/* check to see if we're outside fluid model */
    B = get_model_b(X);		/* field in G */

    if (Ne <= 0.){
      *jnuinv = 0.;
      *knuinv = 0.;
      return;
    }

    /* get covariant four-velocity of fluid for use in get_bk_angle and get_fluid_nu */
    get_model_ucov(X, Ucov);
    get_model_bcov(X, Bcov);
    /*extra: print out stuff to test tetrads */
    get_model_ucon(X, Ucon);
    get_model_bcon(X, Bcon);

    gcov_func(X, gcov);
    lower(Kcon, gcov, Kcov);

    theta = get_bk_angle(X, Kcon, Ucov);	/* angle between k & b */
    if (theta <= 0. || theta >= M_PI) {	/* no emission along field */
	*jnuinv = 0.;
	*knuinv = 0.;
	return;
    }

    nu = get_fluid_nu(Kcon, Ucov);	 /* freq in Hz */
    Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
    
    /* assume emission is thermal */
    Bnuinv = Bnu_inv(nu, Thetae);
    *jnuinv = jnu_inv(nu, Thetae, Ne, B, theta);

    if (Bnuinv < SMALL)
      	*knuinv = SMALL;
    else
      *knuinv = *jnuinv / Bnuinv;

    if (isnan(*jnuinv) || isnan(*knuinv)) {
	fprintf(stderr, "\nisnan get_jkinv\n");
	fprintf(stderr, ">> %g %g %g %g %g %g %g %g\n", *jnuinv, *knuinv,
		Ne, theta, nu, B, Thetae, Bnuinv);
    }

    return;

}

