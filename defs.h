/* some coordinate parameters */
double a;
double R0 ;
double Rin ;
double Rout ;
double Rh,Risco ;
double hslope ;
double th_len, th_beg;
double startx[NDIM], stopx[NDIM], dx[NDIM];
double gam,game,gamp ;
double DTd;
double t0;

double sigma_cut;

//additional for fmks
double mks_smooth,poly_alpha,poly_norm,poly_xt;

/* electron temperature model parameters */
double trat_j;
double theta_j;
double trat_d;
double rmax;

/* HARM model globals */
double M_unit;
double L_unit;
double T_unit;
double RHO_unit;
double U_unit;
double B_unit;

int N1, N2, N3;

double ****bcon;
double ****bcov;
double ****ucon;
double ****ucov;
double ****p;
double ***ne;
double ***uu;
double ***thetae;
double ***b;
double ***sigma_m;


