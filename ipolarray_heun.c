/**********************************************************/
/*** all you need to make a polarized radiative transfer***/
/******* used in ipole to evolve complex tensor N *********/
/***** along with standard evolution for I scalar *********/
/**********************************************************/
/**** written by Monika Moscibrodzka on 09 July 2014 ******/
/************  last update: 6 Oct 2017   ******************/
/****************** co-author: C.F. Gammie ****************/
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

	//Huen's method assuming constant coefficients 
	double SIt,SQt,SUt,SVt;

	// intermediate value 
	SIt = dlam*jI - dlam*(aI*SI0 + aQ*SQ0 + aU*SU0 + aV*SV0) + SI0;
	SQt = dlam*jQ - dlam*(aQ*SI0 + aI*SQ0 + rV*SU0 - rU*SV0) + SQ0;
	SUt = dlam*jU - dlam*(aU*SI0 - rV*SQ0 + aI*SU0 + rQ*SV0) + SU0;
	SVt = dlam*jV - dlam*(aV*SI0 + rU*SQ0 - rQ*SU0 + aI*SV0) + SV0;
		
	// final approximation
	SI = dlam*jI - 0.5*dlam*(aI*SI0 + aQ*SQ0 + aU*SU0 + aV*SV0 + aI*SIt + aQ*SQt + aU*SUt + aV*SVt) + SI0;
	SQ = dlam*jQ - 0.5*dlam*(aQ*SI0 + aI*SQ0 + rV*SU0 - rU*SV0 + aQ*SIt + aI*SQt + rV*SU0 - rU*SVt) + SQ0;
	SU = dlam*jU - 0.5*dlam*(aU*SI0 - rV*SQ0 + aI*SU0 + rQ*SV0 + aU*SIt - rV*SQt + aI*SUt + rQ*SVt) + SU0;
	SV = dlam*jV - 0.5*dlam*(aV*SI0 + rU*SQ0 - rQ*SU0 + aI*SV0 + aV*SIt + rU*SQt - rQ*SUt + aI*SVt) + SV0;
	
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
