#include "../paul.h"

//Relativistic Hydro using TM EoS from Mignone+ 2005

static double RHO_FLOOR = 0.0;
static double PRE_FLOOR = 0.0;

void setHydroParams( struct domain * theDomain ){
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
}

double get_vr( double * prim ){
   double ur = prim[UU1];
   double ut = prim[UU2];
   double u0 = sqrt( 1. + ur*ur + ut*ut );
   return( ur/u0 );
}

double get_h(double Th)
{
    // Th = P / rho
    
    return 2.5*Th + sqrt(1 + 2.25*Th*Th);
}

double get_cs2(double Th, double h)
{
    return Th * (5*h-8*Th) / (3*h*(h-Th));
}

double get_entropy( double * prim )
{
    double Th = prim[PPP] / prim[RHO];
    double h = get_h(Th);
    double es = prim[PPP] * (h-Th) * pow(prim[RHO], -5.0/3.0);
    return log( es); 
}

void prim2cons( double * prim , double * cons , double r , double th , double dV ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double ur  = prim[UU1];
   double ut  = prim[UU2];
   double u0  = sqrt( 1. + ur*ur + ut*ut );

   double rhoh = rho * get_h(prim[PPP] / prim[RHO]);

   cons[DEN] = rho*u0*dV;
   cons[SS1] = rhoh*u0*ur*dV;
   cons[SS2] = r*rhoh*u0*ut*dV;
   cons[TAU] = (rhoh*u0*u0 - Pp - rho*u0)*dV;

   int q;
   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q ){
      cons[q] = cons[DEN]*prim[q];
   }

}

void newt_f(double P, double *f, double *df, double D, double S2, double E)
{
    double v2 = S2 / ((E+P)*(E+P));
    double g = 1.0 / sqrt(fabs(1-v2));
    double Th = g*P / D;
    double h = get_h(Th);

    double dgdP = -v2 * g*g*g / (E+P);
    double dThdP = (g + P*dgdP) / D;
    double dhdP = (2.5*h - 4*Th) / (h - 2.5*Th) * dThdP;

    *f  = D*h*g - E - P;
    *df = D*(dhdP*g + h*dgdP) - 1.0;
}

void cons2prim( double * cons , double * prim , double r , double th , double dV ){

   double D   = cons[DEN]/dV;
   double Sr  = cons[SS1]/dV;
   double St  = cons[SS2]/dV/r;
   double tau = cons[TAU]/dV;
   double E = tau+D;
   double S2 = Sr*Sr + St*St;

   double Pguess = prim[PPP];
   double f,dfdp;
   newt_f( Pguess , &f , &dfdp , D , S2 , E );
   int stop=0;
   while( fabs(f/Pguess) > 1e-8 && stop<10 ){
      Pguess -= f/dfdp;
      newt_f( Pguess , &f , &dfdp , D , S2 , E );
      ++stop;
   }
   
   double Pp = Pguess;
   double v2  = S2/((E + Pp)*(E + Pp));
   double rho = D*sqrt(fabs(1.-v2));
   double vr = Sr/(E+Pp);
   double vt = St/(E+Pp);
   double ur = vr/sqrt(fabs(1.-v2));
   double ut = vt/sqrt(fabs(1.-v2));

   if( rho< RHO_FLOOR ) rho=RHO_FLOOR;
   if( Pp < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[UU1] = ur;
   prim[UU2] = ut;

   int q;
   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q ){
      prim[q] = cons[q]/cons[DEN];
   }

}

void getUstar( double * prim , double * Ustar , double r , double th , double Sk , double Ss , double * n ){

   double rho = prim[RHO];
   double ur  = prim[UU1];
   double ut  = prim[UU2];
   double Pp  = prim[PPP];
   double rhoh = rho * get_h(prim[PPP] / prim[RHO]);

   double un = ur*n[0] + ut*n[1];
   double u0  = sqrt(1.+ur*ur+ut*ut);
   double vn = un/u0;

   double kappa = (Sk-vn)/(Sk-Ss);
   double mn = rhoh*u0*un;
   double E  = rhoh*u0*u0 - Pp;

   double Ps = ( Pp - ( Sk + Ss - vn )*mn + Sk*Ss*E )/( 1.0 - Sk*Ss );

   double alpha1 = ( Ps - Pp )/( Sk - Ss );
   double alpha2 = ( Ss*Ps - vn*Pp )/( Sk - Ss );

   Ustar[DEN] = kappa*rho*u0;
   Ustar[SS1] = kappa*rhoh*u0*ur + alpha1*n[0];
   Ustar[SS2] = r*( kappa*rhoh*u0*ut + alpha1*n[1] );
   Ustar[TAU] = kappa*E + alpha2 - kappa*rho*u0;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DEN];
   }

}

void flux( double * prim , double * flux , double r , double th , double * n ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double ur  = prim[UU1];
   double ut  = prim[UU2];
   double u0  = sqrt(1.+ur*ur+ut*ut);
   double un  = ur*n[0] + ut*n[1];
   double rhoh = rho * get_h(prim[PPP] / prim[RHO]);
 
   flux[DEN] = rho*un;
   flux[SS1] = rhoh*ur*un + Pp*n[0];
   flux[SS2] = r*( rhoh*ut*un + Pp*n[1] );
   flux[TAU] = (rhoh*u0-rho)*un;

   int q;
   for( q=NUM_C ; q<NUM_C+NUM_N ; ++q ){
      flux[q] = flux[DEN]*prim[q];
   }

}

void source(double *prim, double *cons, double *xp, double *xm, double dVdt)
{
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double rhoh = rho * get_h(prim[PPP] / prim[RHO]);
    double rp = xp[0];
    double rm = xm[0];
    double r  = .5*(rp+rm);
    double th = .5*(xp[1]+xm[1]);
    double ut  = prim[UU2];
    double r2 = (rp*rp+rm*rm+rp*rm)/3.;
    cons[SS1] += (2.*Pp + rhoh*ut*ut)*(r/r2)*dVdt;
    cons[SS2] += Pp*(cos(th)/sin(th))*dVdt;
}

void vel( double * prim1 , double * prim2 , double * Sl , double * Sr , double * Ss , double * n , double r , double th ){
   
   double ur1  = prim1[UU1];
   double ut1  = prim1[UU2];
   double un1  = ur1*n[0]+ut1*n[1];
   double W1   = sqrt(1.+ur1*ur1+ut1*ut1);
   double vn1  = un1/W1;

   double rho  = prim1[RHO];
   double Pp1  = prim1[PPP];
   double Th = Pp1 / rho;
   double h = get_h(Th);
   double rhoh1 = rho * h;
   double cs = sqrt(get_cs2(Th, h));
 
   double sigmas = cs*cs/(W1*W1*(1.0-cs*cs));
   *Sl = (vn1 - sqrt( sigmas*(1.0 - vn1*vn1 + sigmas) ) )/( 1.0 + sigmas );
   *Sr = (vn1 + sqrt( sigmas*(1.0 - vn1*vn1 + sigmas) ) )/( 1.0 + sigmas );
 
   double ur2  = prim2[UU1];
   double ut2  = prim2[UU2];
   double un2  = ur2*n[0]+ut2*n[1];
   double W2   = sqrt(1.+ur2*ur2+ut2*ut2);
   double vn2  = un2/W2;

   rho  = prim2[RHO];
   double Pp2  = prim2[PPP];
   Th = Pp2 / rho;
   h = get_h(Th);
   double rhoh2 = rho * h;
   cs = sqrt(get_cs2(Th, h));

   sigmas = cs*cs/(W2*W2*(1.0-cs*cs));
   double sl = (vn2 - sqrt( sigmas*(1.0 - vn2*vn2 + sigmas) ) )/( 1.0 + sigmas );
   double sr = (vn2 + sqrt( sigmas*(1.0 - vn2*vn2 + sigmas) ) )/( 1.0 + sigmas );

   if( *Sl > sl ) *Sl = sl;
   if( *Sr < sr ) *Sr = sr;

   double El = rhoh1*W1*W1 - Pp1;
   double Er = rhoh2*W2*W2 - Pp2;
   double Ml = rhoh1*W1*un1;
   double Mr = rhoh2*W2*un2;
   double Fl = rhoh1*un1*un1 + Pp1;
   double Fr = rhoh2*un2*un2 + Pp2;

   double aL = *Sr;
   double aR = -*Sl;

   double FE = ( aL*Ml + aR*Mr + aL*aR*( El - Er ) )/(aL + aR); 
   double UE = ( aR*El + aL*Er + Ml - Mr )/(aL + aR); 
   double FM = ( aL*Fl + aR*Fr + aL*aR*( Ml - Mr ) )/(aL + aR); 
   double UM = ( aR*Ml + aL*Mr + Fl - Fr )/(aL + aR); 

   if( fabs(FE*UM/pow(UE+FM,2.0)) < 1e-10 ){
      *Ss = UM/(UE + FM); 
   }else{
      *Ss = ( UE + FM - sqrt(fabs((UE + FM)*(UE + FM) - 4.0*FE*UM )) )/( 2.0*FE );
   }
   
}

double get_dL( double * , double * , int );

double get_maxv( double * prim , double w , int dim ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double ur  = prim[UU1];
   double ut  = prim[UU2];
   double W2  = 1. + ur*ur + ut*ut;
   double vr  = ur/sqrt(W2);
   double vt  = ut/sqrt(W2);
   double Th = Pp/rho;
   double h = get_h(Th);

   double cs2 = get_cs2(Th, h);
   double sigmas = cs2/(W2*(1.0-cs2));

   double vn = vr;
   if( dim==1 ) vn = vt;

   double sl = (vn - sqrt( sigmas*(1.0 - vn*vn + sigmas) ) )/( 1.0 + sigmas );
   double sr = (vn - sqrt( sigmas*(1.0 - vn*vn + sigmas) ) )/( 1.0 + sigmas );

   if( dim==0 ){
      sl -= w;
      sr -= w;
   }
   double maxv = fabs(sl);
   if( maxv < fabs(sr) ) maxv = fabs(sr);

   return(maxv);
}

double mindt( double * prim , double w , double * xp , double * xm ){

   double dtr = get_dL(xp,xm,0)/get_maxv(prim,w,0);
   double dtt = get_dL(xp,xm,1)/get_maxv(prim,w,1);;

   double dt = dtr;
   if( dt > dtt ) dt = dtt;

   return( dt );

}
