
#include "../paul.h"

static double T_MIN = 0.0; 
static double adiabatic_index = 0.0;

void setICparams( struct domain * theDomain ){
   T_MIN = theDomain->theParList.t_min;
   adiabatic_index = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double * x ){

   double r  = x[0];
   double th = x[1];

   double t   = T_MIN;
   double th_j = M_PI;

   double Gam = 100.0;

   double E   = 1.0e0;
   double rho0 = 1.0;
   double P0 = 1.0e-5;
   //double r_Sedov = pow(E/rho0,1./3.);


   double betaM = sqrt((1-1/Gam)*(1+1/Gam));

   double R = t * betaM;
   double rhomax = E/(2*M_PI*t*t*t);

   double chi = 1.0 / ((1+betaM)*Gam*Gam * (1-r/t)); // (1-R/t)/(1-r/t)

   double rho, v, X;

   if( r>R || th>th_j ){
      X = 0.0;
      rho = rho0;
      v = 0.0;
   }
   else
   {
       X = 1.0;
       rho = rhomax * chi;
       v = r/t;
   }

   double u = v/sqrt((1-v)*(1+v));
   double Pp = P0 * pow(rho/rho0, adiabatic_index);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = u;
   prim[UU2] = 0.0;

   if( NUM_N > 0 ) prim[NUM_C] = X;

}
