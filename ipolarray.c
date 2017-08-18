
/**********************************************************/

/*** all you need to make a polarized radiative transfer***/
/***** used in ipole to evolve complex tensor N *******/
/***** along with standard evolution for I scalar *********/
/**********************************************************/
/**** written by Monika Moscibrodzka on 09 July 2014 ******/
/************ @ Eindhoven Airport *************************/
/************  last update: 9 May 2017   ******************/

/****************and then rewritten by C.Gammie ***********/

/**********************************************************/

#include "decs.h"
#include "defs.h"


/* the following definitions are used only locally */
#define MNLOOP   for(m=0;m<NDIM;m++)for(n=0;n<NDIM;n++)

/* transfer coefficients in tetrad frame */
void jar_calc(double X[NDIM], double Kcon[NDIM],
	      double *jI, double *jQ, double *jU, double *jV,
	      double *aI, double *aQ, double *aU, double *aV,
	      double *rQ, double *rU, double *rV);

/* tensor tools*/
void check_N(double complex N[NDIM][NDIM], double Kcon[NDIM],
	     double gcov[NDIM][NDIM]);
void complex_lower(double complex N[NDIM][NDIM], double gcov[NDIM][NDIM],
		   int low1, int low2, double complex Nl[NDIM][NDIM]);
void stokes_to_tensor(double fI, double fQ, double fU, double fV,
		      double complex f_tetrad[NDIM][NDIM]);
void tensor_to_stokes(double complex f_tetrad[NDIM][NDIM], double *fI,
		      double *fQ, double *fU, double *fV);
void complex_coord_to_tetrad_rank2(double complex T_coord[NDIM][NDIM],
				   double Ecov[NDIM][NDIM],
				   double complex T_tetrad[NDIM][NDIM]);
void complex_tetrad_to_coord_rank2(double complex T_tetrad[NDIM][NDIM],
				   double Econ[NDIM][NDIM],
				   double complex T_coord[NDIM][NDIM]);

/***************************MAIN FUNCTIONS******************************/
/* initialize tensor N in the coordinate frame at the bening of the *
geodesics integration = it is zero */
void init_N(double X[NDIM], double Kcon[NDIM],
	    double complex N_coord[NDIM][NDIM])
{
    int m, n;

    MNLOOP N_coord[m][n] = 0.0 + I * 0.0;

    return;
}


/*

    parallel transport N over dl 

*/
void push_polar(double Xi[NDIM], double Xm[NDIM], double Xf[NDIM],
		double Ki[NDIM], double Km[NDIM], double Kf[NDIM],
		complex double Ni[NDIM][NDIM],
		complex double Nm[NDIM][NDIM],
		complex double Nf[NDIM][NDIM], double dl)
{

    /* find the connection */
    double lconn[NDIM][NDIM][NDIM];
    get_connection(Xm, lconn);
    int i, j, k, l;

    /* push N */
    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    Nf[i][j] = Ni[i][j];

    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    for (k = 0; k < 4; k++)
		for (l = 0; l < 4; l++)
		    Nf[i][j] += -(lconn[i][k][l] * Nm[k][j] * Km[l] +
				  lconn[j][k][l] * Nm[i][k] * Km[l]
			) * dl / (L_unit * HPL / (ME * CL * CL));

    return;
}

/* updates N for one step on geodesics, using the previous step N*/
/* here we compute new right-hand side of the equation */
/* and somehow rotate this along the geodesics knowing */
/* first point and last point X and K*/

void evolve_N(double Xi[NDIM], double Kconi[NDIM],
	      double Xhalf[NDIM], double Kconhalf[NDIM],
	      double Xf[NDIM], double Kconf[NDIM],
	      double dlam, double complex N_coord[NDIM][NDIM])
{
    int k;
    double gcov[NDIM][NDIM];
    double Ucon[NDIM],Bcon[NDIM];
    double Ecov[NDIM][NDIM], Econ[NDIM][NDIM];
    double complex Nh[NDIM][NDIM];
    double complex N_tetrad[NDIM][NDIM];
    double B;
    double jI, jQ, jU, jV;
    double aI, aQ, aU, aV;
    double rV, rU, rQ; 
    double SI, SQ, SU, SV;
    double SI0, SQ0, SU0, SV0;

    int radiating_region(double X[4]);

    /* parallel transport N by a half, and then full, step */
    push_polar(Xi, Xi, Xhalf, Kconi, Kconi, Kconhalf, N_coord, N_coord, Nh, 0.5 * dlam);
    push_polar(Xi, Xhalf, Xf, Kconi, Kconhalf, Kconf, N_coord, Nh, N_coord, dlam);

    /* absorption/emission/rotation step.  only complete if radiating_region condition is satisfied */
    if ( radiating_region(Xf) ) {

	/* evaluate transport coefficients */
	gcov_func(Xf, gcov);
	jar_calc(Xf, Kconf, &jI, &jQ, &jU, &jV,
		 &aI, &aQ, &aU, &aV, &rQ, &rU, &rV);

	/* make plasma tetrad */
	get_model_ucon(Xf, Ucon);
	B = get_model_b(Xf);	/* field in G */
	if (B > 0.) {
	    get_model_bcon(Xf, Bcon);
	}
	else {
	    Bcon[0] = 0.;
	    for (k = 1; k < NDIM; k++)
		Bcon[k] = 1.;
	}
	make_plasma_tetrad(Ucon, Kconf, Bcon, gcov, Econ, Ecov);

	/* convert N to Stokes */
	complex_coord_to_tetrad_rank2(N_coord, Ecov, N_tetrad);
	tensor_to_stokes(N_tetrad, &SI0, &SQ0, &SU0, &SV0);

	int k,l;
	double M1[NDIM][NDIM];
	double M2[NDIM][NDIM];
	double M3[NDIM][NDIM];
	double M4[NDIM][NDIM];
	double O[NDIM][NDIM];
	double P[NDIM][NDIM];

	double alpha2=aQ*aQ + aU*aU + aV*aV;
	double rho2  =rQ*rQ + rU*rU + rV*rV;
	double alphadrho = aQ*rQ + aU*rU + aV*rV;
	double sig=copysign(1.,alphadrho);
	double T=2.*sqrt(pow(alpha2 - rho2,2)/4.+ alphadrho*alphadrho ) ;
	double ith=1./(T+SMALL); //so that if a=r=0  M234->0, not NaN
	double L1=sqrt(T*0.5+(alpha2-rho2)*0.5)+SMALL; //so that if a=0 alone fac1 != 1/0
	double L2=sqrt(T*0.5-(alpha2-rho2)*0.5)+SMALL; //so that if a=r=0 fac2 != 0   

	for(k=0;k<NDIM;k++)for(l=0;l<NDIM;l++) M1[k][l]=0.0;
	M1[0][0]=1.0;
	M1[1][1]=1.0;
	M1[2][2]=1.0;
	M1[3][3]=1.0;

	M2[0][0]=0.0;
	M2[0][1]=ith*(L2*aQ - sig*L1*rQ);
	M2[0][2]=ith*(L2*aU - sig*L1*rU);
	M2[0][3]=ith*(L2*aV - sig*L1*rV);

	M2[1][0]=ith*(L2*aQ - sig*L1*rQ);
	M2[1][1]=0.0;
	M2[1][2]=ith*(sig*L1*aV + L2*rV);
	M2[1][3]=ith*(-sig*L1*aU - L2*rU);

	M2[2][0]=ith*(L2*aU - sig*L1*rU);
	M2[2][1]=ith*(-sig*L1*aV - L2*rV);
	M2[2][2]=0.0;
	M2[2][3]=ith*(sig*L1*aQ + L2*rQ);

	M2[3][0]=ith*(L2*aV - sig*L1*rV);
	M2[3][1]=ith*(sig*L1*aU + L2*rU);
	M2[3][2]=ith*(-sig*L1*aQ - L2*rQ);
	M2[3][3]=0.0;

	M3[0][0]=0.0;
	M3[0][1]=ith*(L1*aQ + sig*L2*rQ);
	M3[0][2]=ith*(L1*aU + sig*L2*rU);
	M3[0][3]=ith*(L1*aV + sig*L2*rV);

	M3[1][0]=ith*(L1*aQ + sig*L2*rQ);
	M3[1][1]=0.0;
	M3[1][2]=ith*(-sig*L2*aV + L1*rV);
	M3[1][3]=ith*(sig*L2*aU - L1*rU);

	M3[2][0]=ith*(L1*aU + sig*L2*rU);
	M3[2][1]=ith*(sig*L2*aV - L1*rV);
	M3[2][2]=0.0;
	M3[2][3]=ith*(-sig*L2*aQ + L1*rQ);

	M3[3][0]=ith*(L1*aV + sig*L2*rV);
	M3[3][1]=ith*(-sig*L2*aU + L1*rU);
	M3[3][2]=ith*(sig*L2*aQ - L1*rQ);
	M3[3][3]=0.0;


	M4[0][0]=2.*ith*(alpha2+rho2)/2.;
	M4[0][1]=2.*ith*(aV*rU - aU*rV);
	M4[0][2]=2.*ith*(aQ*rV - aV*rQ);
	M4[0][3]=2.*ith*(aU*rQ - aQ*rU);

	M4[1][0]=2.*ith*(aU*rV - aV*rU);
	M4[1][1]=2.*ith*(aQ*aQ + rQ*rQ - (alpha2+rho2)/2.);
	M4[1][2]=2.*ith*(aQ*aU + rQ*rU);
	M4[1][3]=2.*ith*(aV*aQ + rV*rQ);

	M4[2][0]=2.*ith*(aV*rQ - aQ*rV);
	M4[2][1]=2.*ith*(aQ*aU + rQ*rU);
	M4[2][2]=2.*ith*(aU*aU + rU*rU - (alpha2+rho2)/2.);
	M4[2][3]=2.*ith*(aU*aV + rU*rV);

	M4[3][0]=2.*ith*(aQ*rU - aU*rQ);
	M4[3][1]=2.*ith*(aV*aQ + rV*rQ);
	M4[3][2]=2.*ith*(aU*aV + rU*rV);
	M4[3][3]=2.*ith*(aV*aV + rV*rV - (alpha2+rho2)/2.);

	double fac1=1./(aI*aI - L1*L1);
	double fac2=1./(aI*aI + L2*L2);
	double l1dlam=L1*dlam;
	double l2dlam=L2*dlam;
	double EaIdlam=exp(-aI*dlam);
	double coshl1dlam=cosh(l1dlam);
	double sinhl1dlam=sinh(l1dlam);
	double cosl2dlam=cos(l2dlam);
	double sinl2dlam=sin(l2dlam);

	for(k=0;k<4;k++)for(l=0;l<4;l++){

	    O[k][l] = EaIdlam*( 0.5*(coshl1dlam+cosl2dlam)*M1[k][l]
				- sinl2dlam*M2[k][l]
				- sinhl1dlam*M3[k][l]
				+ 0.5*(coshl1dlam-cosl2dlam)*M4[k][l]
				);

	    P[k][l] = (-L1*fac1*M3[k][l] + 0.5*aI*fac1*(M1[k][l] + M4[k][l])) +
	              (-L2*fac2*M2[k][l] + 0.5*aI*fac2*(M1[k][l] - M4[k][l])) -
	      EaIdlam*(
		       ( -L1*fac1*M3[k][l] + 0.5*aI*fac1*(M1[k][l] + M4[k][l]) )*coshl1dlam +
		       ( -L2*fac2*M2[k][l] + 0.5*aI*fac2*(M1[k][l] - M4[k][l]) )*cosl2dlam  +
		       ( -aI*fac2*M2[k][l] - 0.5*L2*fac2*(M1[k][l] - M4[k][l]) )*sinl2dlam  -
		       (  aI*fac1*M3[k][l] - 0.5*L1*fac1*(M1[k][l] + M4[k][l]) )*sinhl1dlam
		      );

	  }

  	 SI = P[0][0]*jI + P[0][1]*jQ + P[0][2]*jU + P[0][3]*jV + O[0][0]*SI0 + O[0][1]*SQ0 + O[0][2]*SU0 + O[0][3]*SV0;
	 SQ = P[1][0]*jI + P[1][1]*jQ + P[1][2]*jU + P[1][3]*jV + O[1][0]*SI0 + O[1][1]*SQ0 + O[1][2]*SU0 + O[1][3]*SV0;
	 SU = P[2][0]*jI + P[2][1]*jQ + P[2][2]*jU + P[2][3]*jV + O[2][0]*SI0 + O[2][1]*SQ0 + O[2][2]*SU0 + O[2][3]*SV0;
	 SV = P[3][0]*jI + P[3][1]*jQ + P[3][2]*jU + P[3][3]*jV + O[3][0]*SI0 + O[3][1]*SQ0 + O[3][2]*SU0 + O[3][3]*SV0;
 
	/* re-pack the Stokes parameters into N */
	stokes_to_tensor(SI, SQ, SU, SV, N_tetrad);
	complex_tetrad_to_coord_rank2(N_tetrad, Econ, N_coord);

    }

    /* SOURCE STEP DONE */

}


/* converts tensor N to Stokes parameters detected at the camera*/
void project_N(double X[NDIM], double Kcon[NDIM],
	       double complex N_coord[NDIM][NDIM], double *Stokes_I,
	       double *Stokes_Q, double *Stokes_U, double *Stokes_V)
{
    double complex N_tetrad[NDIM][NDIM];
    double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];

    make_camera_tetrad(X, Econ, Ecov);

    complex_coord_to_tetrad_rank2(N_coord, Ecov, N_tetrad);

    tensor_to_stokes(N_tetrad, Stokes_I, Stokes_Q, Stokes_U, Stokes_V);

    return;

}

/***************************END MAIN FUNCTIONS******************************/


/*************************SUPPORTING FUNCTIONS******************************/

/*

    call this function if you want to check
    that the coherency tensor N satisfies certain
    basic properties:
    k . N = N . k = 0
    hermitian
    evaluate the invariants: I, Q^2 + U^2, V^2

*/

void check_N(double complex N[NDIM][NDIM],
	     double Kcon[NDIM], double gcov[NDIM][NDIM])
{
    double complex dot;
    double Kcov[NDIM];
    int i, j;

    fprintf(stderr, "enter check_N\n");

    /* compute k . N */
    lower(Kcon, gcov, Kcov);
    fprintf(stderr, "(k . N\n");
    /* first one way */
    for (i = 0; i < 4; i++) {
	dot = 0. + I * 0.;
	for (j = 0; j < 4; j++)
	    dot += Kcov[j] * N[j][i];
	fprintf(stderr, "%d %g + i %g\n", i, creal(dot), cimag(dot));
    }
    /* then the other */
    for (i = 0; i < 4; i++) {
	dot = 0. + I * 0.;
	for (j = 0; j < 4; j++)
	    dot += Kcov[j] * N[i][j];
	fprintf(stderr, "%d %g + i %g\n", i, creal(dot), cimag(dot));
    }
    fprintf(stderr, "k . N)\n");

    /* check for hermiticity */
    fprintf(stderr, "(herm:\n");
    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    fprintf(stderr, "%d %d %g + i %g\n", i, j,
		    creal(N[i][j] - conj(N[j][i])),
		    cimag(N[i][j] - conj(N[j][i]))
		);
    fprintf(stderr, "herm)\n");

    /* check invariants */
    double complex Nud[NDIM][NDIM];
    void complex_lower(double complex N[NDIM][NDIM],
		       double gcov[NDIM][NDIM], int low1, int low2,
		       double complex Nl[NDIM][NDIM]);
    complex_lower(N, gcov, 0, 1, Nud);
    for (i = 0; i < 4; i++)
	fprintf(stderr, "N: %d %g + i %g\n", i, creal(N[i][i]),
		cimag(N[i][i]));
    for (i = 0; i < 4; i++)
	fprintf(stderr, "Nud: %d %g + i %g\n", i, creal(Nud[i][i]),
		cimag(Nud[i][i]));
    dot = 0. + I * 0.;
    for (i = 0; i < 4; i++)
	dot += Nud[i][i];
    fprintf(stderr, "I: %g + i %g\n", creal(dot), cimag(dot));

    double complex Ndd[NDIM][NDIM];
    complex_lower(N, gcov, 1, 1, Ndd);
    dot = 0. + I * 0.;
    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    dot +=
		2. * 0.25 * (N[i][j] + N[j][i]) * (Ndd[i][j] + Ndd[j][i]);
    fprintf(stderr, "IQUsq: %g + i %g\n", creal(dot), cimag(dot));

    dot = 0. + I * 0.;
    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    dot +=
		-2. * 0.25 * (N[i][j] - N[j][i]) * (Ndd[i][j] - Ndd[j][i]);
    fprintf(stderr, "Vsqsq: %g + i %g\n", creal(dot), cimag(dot));

    fprintf(stderr, "leave check_N\n");
}

void complex_lower(double complex N[NDIM][NDIM],
		   double gcov[NDIM][NDIM],
		   int low1, int low2, double complex Nl[NDIM][NDIM]
    )
{
    int i, j, k, l;

    if (!low1 && !low2)
	return;

    if (low1 && low2) {
	for (i = 0; i < 4; i++)
	    for (j = 0; j < 4; j++) {
		Nl[i][j] = 0. + I * 0.;
		for (k = 0; k < 4; k++)
		    for (l = 0; l < 4; l++) {
			Nl[i][j] += N[k][l] * gcov[k][i] * gcov[l][j];
		    }
	    }
	return;
    }

    if (low1) {
	for (i = 0; i < 4; i++)
	    for (j = 0; j < 4; j++) {
		Nl[i][j] = 0. + I * 0.;
		for (k = 0; k < 4; k++)
		    Nl[i][j] += N[k][j] * gcov[k][i];
	    }
	return;
    }

    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++) {
	    Nl[i][j] = 0. + I * 0.;
	    for (l = 0; l < 4; l++)
		Nl[i][j] += N[i][l] * gcov[l][j];
	}
    return;

}

void stokes_to_tensor(double fI, double fQ, double fU, double fV,
		      double complex f_tetrad[NDIM][NDIM])
{
    int i, j;

    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    f_tetrad[i][j] = 0. + I * 0.;

    f_tetrad[1][1] = (fI + fQ + 0. * I);
    f_tetrad[1][2] = (fU - I * fV);
    f_tetrad[2][1] = (fU + I * fV);
    f_tetrad[2][2] = (fI - fQ + 0. * I);

}

void tensor_to_stokes(double complex f_tetrad[NDIM][NDIM],
		      double *fI, double *fQ, double *fU, double *fV)
{

    /*here I divide by two to agree with above */
    *fI = creal(f_tetrad[1][1] + f_tetrad[2][2]) / 2;
    *fQ = creal(f_tetrad[1][1] - f_tetrad[2][2]) / 2;
    *fU = creal(f_tetrad[1][2] + f_tetrad[2][1]) / 2;
    *fV = cimag(f_tetrad[2][1] - f_tetrad[1][2]) / 2;

}

void complex_coord_to_tetrad_rank2(double complex T_coord[NDIM][NDIM],
				   double Ecov[NDIM][NDIM],
				   double complex T_tetrad[NDIM][NDIM])
{
    int i, j, k, l;

    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    T_tetrad[i][j] = 0. + I * 0.;

    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    for (k = 0; k < 4; k++)
		for (l = 0; l < 4; l++)
		    T_tetrad[i][j] +=
			T_coord[k][l] * Ecov[i][k] * Ecov[j][l];

    return;
}

void complex_tetrad_to_coord_rank2(double complex T_tetrad[NDIM][NDIM],
				   double Econ[NDIM][NDIM],
				   double complex T_coord[NDIM][NDIM])
{
    int i, j, k, l;

    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    T_coord[i][j] = 0. + I * 0.;

    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    for (k = 0; k < 4; k++)
		for (l = 0; l < 4; l++)
		    T_coord[i][j] +=
			T_tetrad[k][l] * Econ[k][i] * Econ[l][j];

    return;
}

#undef MNLOOP
