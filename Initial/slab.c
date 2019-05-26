
#include "../paul.h"

static double T_MIN = 0.0; 
static double E0 = 0.0;
static double g0 = 0.0;
static double th_j = 0.0;
static double dR = 0.0;

void setICparams( struct domain * theDomain ){
   T_MIN = theDomain->theParList.t_min;
   E0 = theDomain->theParList.Explosion_Energy;
   g0 = theDomain->theParList.Gam_0;
   th_j = theDomain->theParList.Nozzle_th0;
   dR = theDomain->theParList.Nozzle_r0;
}

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];

   double t   = T_MIN;

   double rho0 = 1.0;
   double cs20 = 1.0e-6;

   double Rf = t*sqrt(1.0 - 1.0/(g0*g0));
   double Rb = Rf-dR;

   double rho, u, Pp, X;

   if( r < Rf && r > Rb && th < th_j)
   {
       double V = 4*M_PI/3.0 * (Rf*Rf+Rf*Rb+Rb*Rb) * dR;
       double M = E0 / (g0 - 1.0);
       rho = M / (g0*V);
       Pp = cs20 * rho;
       u = sqrt(g0*g0-1);
       X = 1.0;
   }
   else
   {
       rho = rho0;
       Pp = cs20 * rho0;
       u = 0.0;
       X = 0.0;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = u;
   prim[UU2] = 0.0;

   if( NUM_N > 0 ) prim[NUM_C] = X;

}
