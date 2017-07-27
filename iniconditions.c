/***************************************************************************
 *   Copyright (C) 2008 by Daniel Banuti   *
 *   daniel.banuti@dlr.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <math.h>

#include "iniconditions.h"
#include "thermodynamics.h"
#include "varutils.h"


extern const double thermo_Gamma;
extern const double thermo_R;


void initialconditions(double P[NOEQNS][NOCELLSX], double U[NOEQNS][NOCELLSX])
{

//  init_refraction_shock_Tpu(P);
//  init_refractionramp_shock_Tpu(P);

  init_shocktube_sod(P);

  // init_shocktube_Tpu(P);
  //init_shocktube_rhopu(P);
  //init_nozzle_Tpu(P);
  // init_wave_Tpu(P);
  //init_ambient_Tpu(P);
}




void init_shocktube_sod(double P[NOEQNS][NOCELLSX])
{

  int i;
  int midfield = NOCELLSX/2;

  for(i=0; i<midfield; i++)
  {
    P[P_rho][i]   = 1.0;
    P[P_p][i]     = 1.0;
    P[P_ux][i]    = 0.0;
  }


  for(i=midfield; i<NOCELLSX; i++)
  {
    P[P_rho][i]   = 0.125;
    P[P_p][i]     = 0.1;
    P[P_ux][i]    = 0.0;
  }

  calc_allP_from_rhopu(P);

}



void init_refraction_shock_Tpu(double P[NOEQNS][NOCELLSX])
{

  int i;
  int midfield = NOCELLSX/2;


  double p2p1, a2a1, M2, u2;
  double gm1 = thermo_Gamma-1.;
  double gp1 = thermo_Gamma+1.;
  
  double T1 = 300.;
  double a1 = sqrt(thermo_Gamma*thermo_R*T1);  
  double p1 = 10000.;
  
    
    p2p1 = 10.; 
    
    M2   = sqrt(1.+(thermo_Gamma+1)/(2.*thermo_Gamma)*(p2p1-1.));
    a2a1 = sqrt(1.+2.*gm1/gp1/gp1*(thermo_Gamma*M2-1./M2/M2-gm1));
    u2   = a2a1*a1*M2;
    
//    printf('%f',u2);
    
    


  for(i=0; i<2; i++)
  {
    P[P_T][0]     = P[P_T][1]     = T1;
    P[P_p][0]     = P[P_p][1]     = p1;
    P[P_ux][0]    = P[P_ux][1]    = 290.;
  }


  for(i=2; i<midfield; i++)
  {
    P[P_T][i]     = 300.0;
    P[P_p][i]     = 1000.0;
    P[P_ux][i]    = 0.0;
  }

  for(i=midfield; i<NOCELLSX; i++)
  {
    P[P_T][i]     = 100.0;
    P[P_p][i]     = 1000.0;
    P[P_ux][i]    = 0.0;
  }

  calc_allP_from_Tpu(P);

}



void init_refractionramp_shock_Tpu(double P[NOEQNS][NOCELLSX])
{

  int ramplength = 10;
  int halframp = ramplength/2;
  int i;
  int midfield = NOCELLSX/2;


//incoming shock
  for(i=0; i<2; i++)
  {
    P[P_T][i]     = 300.0;
    P[P_p][i]     = 10000.0;
    P[P_ux][i]    = 0.0;
  }

//region 1
  for(i=2; i<midfield-halframp; i++)
  {
    P[P_T][i]     = 300.0;
    P[P_p][i]     = 1000.0;
    P[P_ux][i]    = 0.0;
  }
  
  
//transition


//ramping = float(i-(midfield-halframp)
double j=0.;
  for(i=midfield-halframp; i<midfield+halframp; i++)
  {
  	
    P[P_T][i]     = 300.0-j*40.;
    P[P_p][i]     = 1000.0;
    P[P_ux][i]    = 0.0;
    j++;
  }

//region 2
  for(i=midfield; i<NOCELLSX; i++)
  {
    P[P_T][i]     = 100.0;
    P[P_p][i]     = 1000.0;
    P[P_ux][i]    = 0.0;
  }

  calc_allP_from_Tpu(P);

}



void init_ambient_Tpu(double P[NOEQNS][NOCELLSX])
{

  int i;

  for (i=0; i<NOCELLSX;i++)
  {
    P[P_p][i]   = 1000.0;
    P[P_rho][i] = 100.0/thermo_Gamma;
    P[P_ux][i]  = 0.0;
  }

  calc_allP_from_rhopu(P);

}











void init_wave_Tpu(double P[NOEQNS][NOCELLSX])
{
///init sinusoidal wave, 1/4 of domain length

  int i;
  double x, k, pmean, pamp, Tmean, gm1rg;

//   pmean   =   100.0;
//   pamp    =     1.0;
//   rhomean =   100.0;

  gm1rg = (thermo_Gamma-1)/thermo_Gamma;
  int wavelength = ILAST/4;
  k = 2*PI/(double)wavelength;

  for (i=0; i<NOCELLSX;i++)
  {
    P[P_p][i]   = 100.0;
    P[P_rho][i] = 100.0/thermo_Gamma;
    P[P_ux][i]  = 0.0;
  }

  for (i=5; i<wavelength+5;i++)
  {
    x = (double)i-5.0;/*
    P[P_p][i] = pamp*sin(k*x)+pmean;
    P[P_T][i] = Tmean*pow((P[P_p][i]/pmean),gm1rg);
    P[P_ux][i]= sqrt(thermo_Gamma*thermo_R*P[P_T][i])-sqrt(thermo_Gamma*thermo_R*Tmean);*/

    P[P_rho][i] = (100.0 - 1.0*sin(k*x))/thermo_Gamma;
    P[P_ux][i]  = -(1.0*sin(k*x))/100.0;
    P[P_p][i]   = (100.0-1.0*sin(k*x));
  }



  calc_allP_from_rhopu(P);
  //calc_allP_from_Tpu(P);
}


void init_shocktube_Tpu(double P[NOEQNS][NOCELLSX])
{

  int i;
  int midfield = NOCELLSX/2;

  for(i=0; i<midfield; i++)
  {
    P[P_T][i]     = 300.0;
    P[P_p][i]     = 10000.0;
    P[P_ux][i]    = 0.0;
  }


  for(i=midfield; i<NOCELLSX; i++)
  {
    P[P_T][i]     = 300.0;
    P[P_p][i]     = 1000.0;
    P[P_ux][i]    = 0.0;
  }

  calc_allP_from_Tpu(P);

}






void init_shocktube_rhopu(double P[NOEQNS][NOCELLSX])
{

  int i;
  int midfield = NOCELLSX/2;

  for(i=0; i<midfield; i++)
  {
    P[P_rho][i]   = 1.0;
    P[P_p][i]     = 1.0;
    P[P_ux][i]    = 0.0;
  }


  for(i=midfield; i<NOCELLSX; i++)
  {
    P[P_rho][i]   = 0.1;
    P[P_p][i]     = 0.1;
    P[P_ux][i]    = 0.0;
  }

  calc_allP_from_rhopu(P);

}


void init_nozzle(double A[NOCELLSX], double dx[NOCELLSX])
{
///constant dx!
  int i;
  double x;

/*
  double a,b,c,d,L;
  double Amax = 10.0;
  double Amin = 1.0;

  L = ILENGTH*dx[0];

  a = 0.5*(Amax-Amin);
  b = 2*PI/L;
  c = -1.5*PI;
  d = 0.5*(Amin+Amax);*/

// 
//   for (i=0; i<NOCELLSX; i++)
//   {
//     A[i] = a*sin(b*(double)i*dx[1]+c) +d;
//  }


for (i=0; i<NOCELLSX; i++)
  {
    x=0.01*(double)i+0.01;
    A[i] = 1.0/x+x*x*x;
 }

/*
for (i=0; i<NOCELLSX; i++)
  {
    x=0.025*(double)i-2.0;
    A[i] = x*x+2.0;
 }*/

}




void init_nozzle_Tpu(double P[NOEQNS][NOCELLSX])
{

  int i;

  for(i=0; i<1; i++)
  {
    P[P_T][i]     = 2000.0;
    P[P_p][i]     = 1000.0;
    P[P_ux][i]    = 0.0;
  }


  for(i=1; i<NOCELLSX; i++)
  {
    P[P_T][i]     = 300.0;
    P[P_p][i]     = 100.0;
    P[P_ux][i]    = 0.0;
  }



  calc_allP_from_Tpu(P);

}







