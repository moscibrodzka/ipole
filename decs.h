#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include "constants.h"
#include <complex.h> 

#ifdef _OPENMP
#include <omp.h>
#endif

/* PARAMETERS TO CHANGE */

//change both to 1 when running Novikov-Thorne model
#define NT_PROB 0
#define SOURCE_NT 0

//or which source SGRA,M87,DABHB, if NT model these all 0
#define SOURCE_SGRA 1
#define SOURCE_M87  0
#define SOURCE_DABHB 0

//image resolution
#define NX  80
#define NY  80

// FOV in [M] 
#define DX 40.
#define DY 40.

//parameters of eDF
#define THERMAL 1
#define POWERL  0
#define KAPPAL  0

//floor and ceiling of electron temperature, model dependent
#define THETAE_MIN 0.001
#define THETAE_MAX 100

#define NPRIM	8 // default
//#define NPRIM	10 // 10 variables in hdf5 in fmks

//chose integration scheme in ipolarray.c
#define INT_FULL  0 //full integration step
#define INT_SPLIT 1 //split: 1/2rot + ea + 1/2rot (from development version of ipole)

//colortheme for ppm files
#define RAINBOW 0
#define AFMHOT 1
#define BW 0


/********************** RATHER DO NOT CHANGE BELOW *************/

#define MAXNSTEP        60000
#define POLARIZATION_ON (1)


/**************************** WARNING *************************/
/********************** DO NOT CHANGE BELOW *******************/

#define NDIM  4

/* mnemonics for primitive vars; conserved vars */
#define KRHO    0
#define UU      1
#define U1      2
#define U2      3
#define U3      4
#define B1      5
#define B2      6
#define B3      7

#define KEL      8
#define KTOT     9

/* numerical convenience */
#define SMALL	1.e-40
#define SMALLER	1.e-50
					     
#define sign(x) (((x) < 0) ? -1 : ((x) > 0))
#define MAX(c,b) (((c)>(b))?(c):(b))
#define MIN(c,b) (((c)<(b))?(c):(b))

/* some coordinate parameters */
extern double a;

extern double theta_j;
extern double trat_j;
extern double trat_d;
extern double rmax;

extern double R0 ;
extern double Rin ;
extern double Rout ;
extern double Rh,Risco ;
extern double hslope ;
extern double th_len,th_beg;
extern double startx[NDIM], stopx[NDIM], dx[NDIM];
extern double gam,game,gamp ;
extern double DTd;
extern double t0;

//additional parameter 
extern double sigma_cut;

//additional global parameters for fmks
extern double mks_smooth,poly_alpha,poly_norm,poly_xt;

/* HARM model globals */
extern double M_unit;
extern double L_unit;
extern double T_unit;
extern double RHO_unit;
extern double U_unit;
extern double B_unit;

extern int N1, N2, N3;

extern double ****bcon;
extern double ****bcov;
extern double ****ucon;
extern double ****ucov;
extern double ****p;
extern double ***ne;
extern double ***uu;
extern double ***thetae;
extern double ***b;



double rhorizon_calc(int pos_sign);
double risco_calc( int do_prograde );

void linear_interp_I_Q_disk(double mu, double *f1, double *f2);
void get_NT_S(double X[NDIM], double Kcon[NDIM], double F, double complex N_coord[NDIM][NDIM]);
double get_NT_nu(double X[NDIM], double Kcon[NDIM]);
double bnu(double nu, double T);
double krolikc(double r, double a, double r_isco);
double get_weight(double *xx, double x, int *jlo);

/** model-independent subroutines **/
/* core routines in main.c */
void init_XK(int i, int j, double Xcam[4], double fovx, double fovy, double X[4], double Kcon[4]) ;
void null_normalize(double Kcon[NDIM], double fnorm) ;
void normalize(double *vcon, double gcov[][NDIM]);
double approximate_solve(double Ii, double ji, double ki, double jf, double kf, double dl) ;
void get_jkinv(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv,double col_v[NDIM]);
int  stop_backward_integration(double X[NDIM], double Kcon[NDIM],double Xcam[NDIM]) ;
void dump(double image[NX][NY],double imageS[NX][NY][NDIM], double col_var[NX][NY][3],char *fname, double scale) ;
void make_txt(double image[NX][NY], double imageS[NX][NY][NDIM], double freqcgs, char *fname, double FOVX, double FOVY,double scale);
double root_find2(double x[NDIM]);


/* geodesic integration */
void push_photon(double X[NDIM], double Kcon[NDIM], double dl, double Xhalf[NDIM],double Kconhalf[NDIM]);

/* model */
void   gcov_func(double *X, double gcov[][NDIM]);
void   gcov_func_rec(double *X, double gcov[][NDIM]);
int   gcon_func(double gcov[][NDIM], double gcon[][NDIM]);
double gdet_func(double gcov[][NDIM]);

void   set_dxdX(double *X, double dxdX[NDIM][NDIM]);
void   get_connection(double *X, double lconn[][NDIM][NDIM]);
void   get_connection_num(double *X, double lconn[][NDIM][NDIM]);

void   lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov);
double stepsize(double X[NDIM], double K[NDIM]);

void init_model(char *args[]) ;
double get_model_thetae(double X[NDIM]) ;
double get_model_b(double X[NDIM]) ;
double get_model_ne(double X[NDIM]) ;
double get_model_uu(double X[NDIM]) ;
void get_model_bcov(double X[NDIM], double Bcov[NDIM]) ;
void get_model_bcon(double X[NDIM], double Bcon[NDIM]) ;
void get_model_ucov(double X[NDIM], double Ucov[NDIM]) ;
void get_model_ucon(double X[NDIM], double Ucon[NDIM]) ;

/* harm utilities */
/*
void interp_fourv(double X[NDIM], double ***fourv, double Fourv[NDIM]) ;
double interp_scalar(double X[NDIM], double **var) ;
void Xtoij(double X[NDIM], int *i, int *j, double del[NDIM]) ;
*/

void  bl_coord(double *X, double *r, double *th);
void set_units(char *instr) ;
void init_physical_quantities(void) ;
void init_storage(void) ;

/* tetrad related */
void   coordinate_to_tetrad(double Ecov[NDIM][NDIM], double K[NDIM], double K_tetrad[NDIM]);
void   tetrad_to_coordinate(double Ecov[NDIM][NDIM], double K_tetrad[NDIM], double K[NDIM]);
double delta(int i, int j);
void   make_camera_tetrad(double X[NDIM], double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM]) ;
void   make_plasma_tetrad(double Ucon[NDIM], double Kcon[NDIM], double Bcon[NDIM],
	double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM]) ;

/* imaging */
void make_ppm(double p[NX][NY], double freq, char filename[]) ;
void afmhot_palette(double data, double min, double max, int *pRed, int *pGreen, int *pBlue) ;
void rainbow_palette(double data, double min, double max, int *pRed, int *pGreen, int *pBlue) ;
void monika_palette(double data, double min, double max, int *pRed, int *pGreen, int *pBlue) ;

/* radiation */
double Bnu_inv(double nu, double Thetae) ;
double jnu_inv(double nu, double Thetae, double Ne, double B, double theta) ;
double get_fluid_nu(double Kcon[NDIM], double Ucov[NDIM]) ;
double get_bk_angle(double X[NDIM], double Kcon[NDIM], double Ucov[NDIM]) ;
double get_bk_angle2(double X[NDIM], double Kcon[NDIM], double Ucov[NDIM], double Bcon[NDIM],double Bcov[NDIM]) ;

/* emissivity */ 
double jnu_synch(double nu, double Ne, double Thetae, double B, double theta) ;
double jnu_synchSI(double nu, double Ne, double Thetae, double B,
		   double theta);
double jnu_synchSQ(double nu, double Ne, double Thetae, double B,
		   double theta);
double jnu_synchSV(double nu, double Ne, double Thetae, double B,
		   double theta);


void init_N(double Xi[NDIM],double Kconi[NDIM],double complex Ncon[NDIM][NDIM]);
void evolve_N(double Xi[NDIM],double Kconi[NDIM],
	      double Xf[NDIM],double Kconf[NDIM],
	      double Xhalf[NDIM],double Kconhalf[NDIM],
	      double dlam,
	      double complex N_coord[NDIM][NDIM],
	      double *tauF);
void project_N(double X[NDIM],double Kcon[NDIM],
	double complex Ncon[NDIM][NDIM],
	double *Stokes_I, double *Stokes_Q,double *Stokes_U,double *Stokes_V);






