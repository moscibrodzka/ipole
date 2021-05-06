#include "decs.h"
#include "defs.h"
#include <omp.h>
#include <time.h>

/*shared*/
double image[NX][NY];
double imageS[NX][NY][NDIM];
double col_var[NX][NY][3];

int main(int argc, char *argv[])
{
    FILE *fp;
    double fovx, fovy;
    double thetacam, phicam, rcam, Xcam[NDIM];
    double freq, freqcgs;
    double Ftot, Dsource;
    double scale;
    int imax, jmax;
    double Imax, Iavg;	     
    double col_v_int[4];

#if(NT_PROB)
    if (argc < 3) {
	fprintf(stderr,"usage: ipole theta freq a\n");
	exit(0);
    }
    sscanf(argv[1], "%lf", &thetacam);
    sscanf(argv[2], "%lf", &freqcgs);
    sscanf(argv[3],"%lf",&a) ;
#else    
    if (argc < 7) {
	fprintf(stderr,"usage: ipole theta freq filename Munit trat_j trat_d sigma_cut\n");
	exit(0);
    }

    sscanf(argv[1], "%lf", &thetacam);
    sscanf(argv[2], "%lf", &freqcgs);
    sscanf(argv[4], "%lf", &M_unit);
    sscanf(argv[5], "%lf", &trat_j);
    sscanf(argv[6], "%lf", &trat_d);
    sscanf(argv[7], "%lf", &sigma_cut);
#endif

    //initiate model and photon frequency
    init_model(argv);
    freq = freqcgs * HPL / (ME * CL * CL);

    //camera position
    rcam = 10000.; // default value
    phicam = 0.0;
    Xcam[0] = 0.0;
    Xcam[1] = log(rcam);
    double x[NDIM] = {0., rcam, thetacam/180.*M_PI, phicam/180.*M_PI};
    Xcam[2] = root_find2(x);
    Xcam[3] = phicam/180.*M_PI;
    fprintf(stdout,"camera coordinates: %g %g %g %g \n",Xcam[0],Xcam[1],Xcam[2],Xcam[3]);
    double r,th;
    bl_coord(Xcam,&r,&th);
    fprintf(stdout,"camera theta=: %g \n",th/M_PI*180.);                                                                                                           
    
#if(SOURCE_SGRA)
    Dsource = DSGRA;
#endif
#if(SOURCE_M87)
    Dsource = DM87;
#endif
#if(SOURCE_DABHB)
    Dsource = DABHB;
#endif
#if(SOURCE_NT)
    Dsource = 0.05*PC;
#endif

    //to chose image FOV in uas one has to uncomment this and problably next
    fovx = DX / rcam;
    fovy = DY / rcam;
    scale = (DX * L_unit / NX) * (DY * L_unit / NY) / (Dsource * Dsource) / JY;
    fprintf(stderr,"a=%g rh=%g \n",a,Rh);
    fprintf(stderr,"intensity [cgs] to flux per pixel [Jy] conversion: %g\n",scale);
    fprintf(stderr,"Dsource: %g [cm]\n",Dsource);
    fprintf(stderr,"Dsource: %g [kpc]\n",Dsource/(1.e3*PC));
    fprintf(stderr,"FOVx, FOVy: %g %g [GM/c^2]\n",DX,DY);
    fprintf(stderr,"FOVx, FOVy: %g %g [rad]\n",DX*L_unit/Dsource,DY*L_unit/Dsource);
    fprintf(stderr,"FOVx, FOVy: %g %g [muas]\n",DX*L_unit/Dsource * 2.06265e11 ,DY*L_unit/Dsource * 2.06265e11);
    
    clock_t begin = clock();


#pragma omp parallel for schedule(dynamic,1) collapse(2) shared(image,imageS,Xcam,fovx,fovy)
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NY; j++) {

	double X[NDIM], Kcon[NDIM], Xhalf[NDIM], Kconhalf[NDIM];    
	struct of_traj {
	  double dl;
	  double X[NDIM];
	  double Kcon[NDIM];
	  double Xhalf[NDIM];
	  double Kconhalf[NDIM];
	} traj[MAXNSTEP];
		
	init_XK(i, j, Xcam, fovx, fovy, X, Kcon);
	for (int k = 0; k < NDIM; k++) Kcon[k] *= freq;
	
	/* integrate geodesic backwards along trajectory */
	int nstep = 0;

	double tauFn,tauI,tauQ,l;
	tauFn=0.0;
	tauI=0.0;
	tauQ=0.0;
	l=0.0;
	double jI, jQ, jU, jV;
	double aI, aQ, aU, aV;
	double rV, rU, rQ;
	int n_ring=0;
	double dxdX[NDIM][NDIM];
	
	while (!stop_backward_integration(X, Kcon, Xcam)){
	  
	  double dl = stepsize(X, Kcon);
	  double X2_p=X[2];
	  push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
	  nstep++;

	  traj[nstep].dl = dl * L_unit * HPL / (ME * CL * CL);
	  for (int l = 0; l < NDIM; l++) traj[nstep].X[l] = X[l];
	  for (int l = 0; l < NDIM; l++) traj[nstep].Kcon[l] = Kcon[l];
	  for (int l = 0; l < NDIM; l++) traj[nstep].Xhalf[l] = Xhalf[l];
	  for (int l = 0; l < NDIM; l++) traj[nstep].Kconhalf[l] = Kconhalf[l];
	  
	  if (nstep > MAXNSTEP - 2) {
	    fprintf(stderr, "MAXNSTEP exceeded on j=%d i=%d\n", j,i);
	    exit(1);
	  }
	  
	}

	nstep--; /* final step violated the "stop" condition,so don't record it */
	/* DONE geodesic integration */

	double Xi[NDIM], Xf[NDIM], Kconi[NDIM], Kconf[NDIM], ji, ki, jf, kf;
	double col_v[NDIM];
	double complex N_coord[NDIM][NDIM];
	double Intensity = 0.0;
	double tauF=0.0;
	
	for (int l = 0; l < NDIM; l++) {
	  Xi[l] = traj[nstep].X[l];
	  Kconi[l] = traj[nstep].Kcon[l];
	}
	init_N(Xi, Kconi, N_coord);

	while (nstep > 1) {
	  for (int l = 0; l < NDIM; l++) {
	    Xi[l]       = traj[nstep].X[l];
	    Kconi[l]    = traj[nstep].Kcon[l];
	    Xhalf[l]    = traj[nstep].Xhalf[l];
	    Kconhalf[l] = traj[nstep].Kconhalf[l];
	    Xf[l]       = traj[nstep - 1].X[l];
	    Kconf[l]    = traj[nstep - 1].Kcon[l];
	  }
	  double dl = traj[nstep].dl; //this is ok ! unless i integrate here again

#if(NT_PROB)
	  //for novicov problem we solve just transport
	  //NT:if geodesics crosses the disk then, assign non-zero I
	  double lrd_in=log(risco_calc(1));
	  double lrd_out=log(100.);
	  if( ( (Xi[2]>0.5 && Xf[2]<0.5) || (Xi[2]<0.5 && Xf[2]>0.5)) && (Xf[1] > lrd_in) && (Xf[1] < lrd_out) ){
	      double r=exp(Xf[1]);
	      double Teff;
	      double Mdotedd = 4. * M_PI * GNEWT * 10.* MSUN* MP / 0.1 / CL / SIGMA_THOMSON;
	      double Mdot = 0.01*Mdotedd;
	      double b = 1. - 3. / r + 2. * a / pow(r, 3. / 2.);
	      double kc = krolikc(r, a, risco_calc(1));
	      double T0=pow(3.0 / 8.0 / M_PI * GNEWT * 10. *MSUN * Mdot / pow(L_unit, 3) / SIG, 1. / 4.);
	      if (r < exp(lrd_out)) {
		  Teff = T0 * pow(kc / b / pow(r, 3), 1. / 4.);
	      } else {
		  Teff = T0 / 1e5;
	      }
	      double ff=1.8;
	      double nu=get_NT_nu(Xf,Kconf); //return nu in fluid frame in Hz, and polarization tensor
	      Intensity=pow(ff,-4)*bnu(nu,ff*Teff);///nu/nu/nu;//invariant Intensity
	      get_NT_S(Xf,Kconf,Intensity,N_coord); //assign correct N_coord in this second round
	      Intensity=Intensity/nu/nu/nu;
	  }
#else	      
	      /* solve total intensity equation alone */
	      get_jkinv(Xi, Kconi, &ji, &ki, col_v);
	      get_jkinv(Xf, Kconf, &jf, &kf, col_v);
	      Intensity = approximate_solve(Intensity, ji, ki, jf, kf, dl);
#endif	  
	
#if(POLARIZATION_ON)
	      /* solve polarized transport */
	      evolve_N(Xi, Kconi, Xhalf, Kconhalf, Xf, Kconf, traj[nstep].dl, N_coord, &tauF);
#endif
	  /* swap start and finish */
	  ji = jf;
	  ki = kf;
	  
	  nstep--;
	}
	  
	image[i][j] = Intensity * pow(freqcgs, 3);
#if(POLARIZATION_ON)
	double Stokes_I,Stokes_Q,Stokes_U,Stokes_V;
	project_N(Xf, Kconf, N_coord, &Stokes_I, &Stokes_Q, &Stokes_U, &Stokes_V);
	imageS[i][j][0] = Stokes_I * pow(freqcgs, 3);
	imageS[i][j][1] = Stokes_Q * pow(freqcgs, 3);
	imageS[i][j][2] = Stokes_U * pow(freqcgs, 3);
	imageS[i][j][3] = Stokes_V * pow(freqcgs, 3);
	col_var[i][j][0]=tauF;
#endif


	if (j==0) fprintf(stderr,"%d ",i); 
	
      } 
    }

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    fprintf(stdout,"\n time spent = %lf \n",time_spent);

    /* printing out to files and on stderr */
    double FtotI=0.0;
    double FtotQ=0.0;
    double FtotU=0.0;
    double FtotV=0.0;
    Ftot = 0.;
    Imax = 0.0;
    Iavg = 0.0;
    imax = jmax = 0;
    for (int i = 0; i < NX; i++)
	for (int j = 0; j < NY; j++) {
	    Ftot += image[i][j] * scale;	

	    FtotI += imageS[i][j][0] * scale;	
	    FtotQ += imageS[i][j][1] * scale;	
	    FtotU += imageS[i][j][2] * scale;	
	    FtotV += imageS[i][j][3] * scale;	

	    Iavg += image[i][j];
	    if (image[i][j] > Imax) {
		imax = i;
		jmax = j;
		Imax = image[i][j];
	    }
	}
    fprintf(stderr, "imax=%d jmax=%d Imax=%g Iavg=%g\n", imax, jmax, Imax,Iavg / (NX * NY));
    fprintf(stderr, "LP =  %g [per cent], CP = %g [per cent]\n",sqrt(FtotQ*FtotQ+FtotU*FtotU)/FtotI*100.,FtotV/FtotI*100.);
    fprintf(stderr, "Ftot: %g %g %g %g scale=%g\n", freqcgs, Ftot,  Ftot * Dsource * Dsource * JY * freqcgs * 4 * M_PI, 
	    Dsource * Dsource * JY  * 4 * M_PI, scale);

    /* image, dump result */
    make_ppm(image, freq, "output_ipole/ipole_fnu.ppm");
    dump(image, imageS, col_var, "output_ipole/ipole.dat", scale);

    //make_txt(image,imageS,freqcgs,"output_ipole/ipole_ehtim.txt",DX*L_unit/Dsource * 206264806247.1,DY*L_unit/Dsource * 206264806247.1,scale);
    for (int i = 0; i < NX; i++)
	for (int j = 0; j < NY; j++)
	    image[i][j] = log(image[i][j] + 1.e-50);
    make_ppm(image, freq, "output_ipole/ipole_lfnu.ppm");
    
    /* done! */

    return 0;
}


/* dump-file for ehtim library */
/* here FOV is in miuarcseconds */
void make_txt(double image[NX][NY],double imageS[NX][NY][NDIM], double freqcgs, char *fname, double FOVX, double FOVY, double scale)
{

  int i, j;
  FILE *fp;

  fp = fopen(fname, "w");
  if (fp == NULL) {
    fprintf(stderr, "Failes to open txt file %s\n", fname);
    exit(125);
  }

  /* write out header information, has to be changed for other sources */
  fprintf(fp,"# SRC: M87\n");
  fprintf(fp,"# RA: 12 h 30 m 49.4234 s\n");
  fprintf(fp,"# DEC: 12 deg 23 m 28.0437 s\n");
  fprintf(fp,"# MJD: 48277\n");
  fprintf(fp,"# RF: %g GHz\n",freqcgs*1e-9);
  fprintf(fp,"# FOVX: %d pix %g as\n",NX,FOVX*1e-6);
  fprintf(fp,"# FOVY: %d pix %g as\n",NY,FOVY*1e-6);
  fprintf(fp,"# ------------------------------------\n");

#if(POLARIZATION_ON)
  fprintf(fp,"# x (as)         y (as)         I (Jy/pixel)     Q (Jy/pixel)     U (Jy/pixel)     \n");
#else
  fprintf(fp,"# x (as)         y (as)         I (Jy/pixel)\n");
#endif
  fflush(fp);
#if(POLARIZATION_ON)
  for (j = NY; j > 0; j--)
    for (i = 0; i < NX; i++)
      {
	fprintf(fp,"%.10f %.10f %g %g %g %g\n",
		(i-NX/2)*FOVX*1e-6/(double)NX,
                (-j+NY/2)*FOVY*1e-6/(double)NY,
                imageS[i][j][0]*scale,
                imageS[i][j][1]*scale,
                imageS[i][j][2]*scale,
		imageS[i][j][3]*scale);
      }
#else
  for (j = NY; j > 0; j--)
    for (i = 0; i < NX; i++)
      {
        fprintf(fp,"%15.10g %15.10g %15.10g\n",
                (i-NX/2)*FOVX*1e-6/(double)NX,
		(-j+NY/2)*FOVY*1e-6/(double)NY,
                image[i][j]*scale);
      }
#endif

  fclose(fp);
  return;
}



void dump(double image[NX][NY], double imageS[NX][NY][NDIM], double col_varp[NX][NY][3],char *fname,
	  double scale)
{
    FILE *fp;
    int i, j;
    double xx, yy, dum_r;
    double sum_i;


    fp = fopen(fname, "w");
    if (fp == NULL) {
	fprintf(stderr, "unable to open %s\n", fname);
	exit(1);
    }


    sum_i = 0.0;
    for (i = 0; i < NX; i++) {
	for (j = 0; j < NY; j++) {
	    sum_i += image[i][j];
	    
	    //rotated QU according to IAU standards with +Q aligned with north
	    fprintf(fp, "%d %d %g %g %g %g %g  %g\n",
		    i, j, 
		    image[i][j]*scale,
		    imageS[i][j][0]*scale, 
		    -imageS[i][j][1]*scale,
		    -imageS[i][j][2]*scale, 
		    imageS[i][j][3]*scale,
		    col_varp[i][j][0]
		    );
	}
    }
    fclose(fp);


    /* dump a vtk file with Stokes I for visit */
    fp = fopen("output_ipole/ipole.vtk", "w");
    if (fp == NULL) {
	fprintf(stderr, "unable to open %s\n", "ipole.vtk");
	exit(1);
    }

    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Image Simulation Result\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", NX , NY , 1);
    fprintf(fp, "POINTS %d float\n", (NX ) * (NY));

    for (j = 0; j < NY; j++) {
	for (i = 0; i < NX; i++) {
	    xx = i * 1.0;	//-x_size/2.+i*stepx;                                                                                                                 
	    yy = j * 1.0;	//-y_size/2.+j*stepy;                                                                                                                 
	    fprintf(fp, "%g %g %g\n", xx, yy, 0.0);
	}
    }
    fprintf(fp, "\nPOINT_DATA %d\n", (NX) * (NY));
    fprintf(fp, "SCALARS Intensity float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    for (j = 0; j < NY; j++) {
	for (i = 0; i < NX; i++) {
	    dum_r = (image[i][j]);
	    fprintf(fp, "%g\n", dum_r);
	}
    }

    fclose(fp);



}

void init_XK(int i, int j, double Xcam[4], double fovx, double fovy,	/* field of view, in radians */
	  double X[4], double Kcon[4]	/* position, wavevector */
    )
{
    double Econ[NDIM][NDIM];
    double Ecov[NDIM][NDIM];
    double Kcon_tetrad[NDIM];
    int k;

    /* construct orthonormal tetrad.
       e^0 along Ucam
       e^3 outward (!) along radius vector
       e^2 toward north pole of coordinate system
       ("y" for the image plane)
       e^1 in the remaining direction
       ("x" for the image plane)

    */
    make_camera_tetrad(Xcam, Econ, Ecov);

    /* construct *outgoing* wavevectors */

    /* EHT library */
    double rotcam=0.0;
    double xoff=0.0;
    double yoff=0.0;

    double dxoff=(xoff + i + 0.5 - 0.01)/NX -0.5;
    double dyoff=(yoff + j + 0.5)/NY -0.5;
    
    Kcon_tetrad[0] = 0.;
    Kcon_tetrad[1] = (dxoff*cos(rotcam)-dyoff*sin(rotcam))*fovx;
    Kcon_tetrad[2] = (dxoff*sin(rotcam)+dyoff*cos(rotcam))*fovy;
    Kcon_tetrad[3] = 1.;
    
    /* normalize */
    null_normalize(Kcon_tetrad, 1.);

    /* translate into coordinate frame */
    tetrad_to_coordinate(Econ, Kcon_tetrad, Kcon);

    /* set position */
    for (k = 0; k < NDIM; k++)
	X[k] = Xcam[k];

    /* done! */
}

/* normalize null vector in a tetrad frame */
void null_normalize(double Kcon[NDIM], double fnorm)
{
    double inorm;

    inorm =
	sqrt(Kcon[1] * Kcon[1] + Kcon[2] * Kcon[2] + Kcon[3] * Kcon[3]);

    Kcon[0] = fnorm;
    Kcon[1] *= fnorm / inorm;
    Kcon[2] *= fnorm / inorm;
    Kcon[3] *= fnorm / inorm;

}

/* 

   must be a stable, approximate solution to radiative transfer
   that runs between points w/ initial intensity I, emissivity
   ji, opacity ki, and ends with emissivity jf, opacity kf.

   Return final intensity

*/

double approximate_solve(double Ii, double ji,
			 double ki, double jf, double kf, double dl)
{
    double efac, If, javg, kavg, dtau;

    javg = (ji + jf) / 2.;
    kavg = (ki + kf) / 2.;

    dtau = dl * kavg;

    if (dtau < 1.e-3) {
	If = Ii + (javg - Ii * kavg) * dl * (1. -
					     (dtau / 2.) * (1. -
							    dtau / 3.));
    } else {
	efac = exp(-dtau);
	If = Ii * efac + (javg / kavg) * (1. - efac);
    }

    
    
    return (If);
}

