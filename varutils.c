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

#include "varutils.h"
#include "thermodynamics.h"

#define MAX(a,b)      ( (a)>=(b) ? (a) : (b) )
#define EPS 1.0e-12



extern const double thermo_Gamma;
extern const double thermo_R;



void reset_flux(double flux[NOEQNS][NOIFX])
{
  int i,eqn;

  for (eqn=0; eqn<NOEQNS; eqn++)
  {
    for (i=0; i<NOIFX; i++)
    {
       flux[eqn][i] = 0.0;
    }
  }
}


void reset_residual(double residual[NOEQNS][NOCELLSX])
{
  int eqn,i;

  for (eqn=0; eqn<NOEQNS; eqn++)
  {
    for (i=0; i<NOCELLSX; i++)
    {
       residual[eqn][i] = 0.0;
    }
  }
}





void initvariables(double        P[NOPRIMS][NOCELLSX],
                   double        U[NOEQNS][NOCELLSX],
                   double residual[NOEQNS][NOCELLSX],
                   double     flux[NOEQNS][NOIFX],
                   double  l_state[NOPRIMS][NOIFX],
                   double  r_state[NOPRIMS][NOIFX],
                   double       dx[NOCELLSX],
                   double        A[NOCELLSX])
{
  int i,j;
  for (i=0; i<NOCELLSX; i++)
  {
    dx[i] = 0.01;
    A[i]  = 1.0;


    for (j=0; j<NOEQNS; j++)
    {
      U[j][i]    = 0.0;
    }

    for (j=0; j<NOPRIMS; j++)
    {
      P[j][i] = 0.0;
      l_state[j][i] = 0.0;
      r_state[j][i] = 0.0;
    }
  }
  reset_flux(flux);
  reset_residual(residual);
}


void cons_prims(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX])
{
  int i;
  for (i=0; i < NOCELLSX; i++)
  {
    P[P_rho][i] = U[U_rho][i];
    P[P_ux][i]  = U[U_rhoux][i] / U[U_rho][i];
    P[P_et][i]  = U[U_rhoet][i] / U[U_rho][i];/*
printf("Cell %i \n",i);
printf("P_rho %f\n",P[P_rho][i]);
printf("P_ux %f\n",P[P_ux][i]);
printf("P_et %f\n",P[P_ux][i]);*/

  }
// printf("Cell %i\n",i);
}



void prims_cons(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX])
{
  int i;
  for (i=0; i < NOCELLSX; i++)
  {
     U[U_rho][i]   = P[P_rho][i];
     U[U_rhoux][i] = P[P_rho][i]*P[P_ux][i];
     U[U_rhoet][i] = P[P_rho][i]*P[P_et][i];
  }
}


double get_timestep(double CFL, double P[NOPRIMS][NOCELLSX], double dx[NOCELLSX])
{

  double dt, umax=0.0, sos, umaxrdx;
  umax,umaxrdx = 0.0;

// printf("\n Hallo Welt in!\n\n");

  for (int i=0; i < NOCELLSX; i++)
  {
  
    sos = sqrt(thermo_Gamma*P[P_p][i]/P[P_rho][i]);
    umaxrdx =  MAX( (sos+fabs(P[P_ux][i]))/dx[i],(umaxrdx));
    // umax =  MAX( (sos+fabs(P[P_ux][i])),(umax));

    // umaxrdx = umax/dx[i];
  }

  // printf("\n Hallo Welt after loop!\n\n %f %f",CFL,umaxrdx);  

  dt = CFL/umaxrdx*0.1;

  // dt = 0.0003;
  // printf('dt %f \n',dt);

  return dt;
}




void fixnegativeprims(double P[NOPRIMS][NOCELLSX])
{
  int i;
  for (i=0; i<NOCELLSX; i++)
  {
  	
  	P[P_rho][i] = MAX(P[P_rho][i],EPS);
  	P[P_p][i]   = MAX(P[P_p][i],EPS);
  	P[P_T][i]   = MAX(P[P_T][i],EPS);
  	
  }
}








