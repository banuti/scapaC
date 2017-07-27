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


#include "thermodynamics.h"
#include "varutils.h"


#define SQR(a) ((a)*(a))


const double thermo_Gamma =   1.4;
const double thermo_R     = 287.0;





void calc_allP_from_Tpu(double P[NOPRIMS][NOCELLSX])
{
  int i;
  for (i=0; i< NOCELLSX; i++)
  {
    P[P_rho][i] = P[P_p][i]/(thermo_R * P[P_T][i]);
    P[P_e][i]   = thermo_R/(thermo_Gamma-1.0)*P[P_T][i];
    P[P_et][i]  = P[P_e][i]+0.5*SQR( P[P_ux][i] );
    P[P_ht][i] = P[P_et][i] + P[P_p][i]/P[P_rho][i];
    P[P_h][i]  = P[P_e][i] + P[P_p][i]/P[P_rho][i];
  }
}



void calc_allstates_from_Tpu(double state[NOPRIMS][NOIFX])
{
  int i;
  for (i=0; i< NOIFX; i++)
  {
    state[P_rho][i] = state[P_p][i]/(thermo_R * state[P_T][i]);
    state[P_e][i]   = thermo_R/(thermo_Gamma-1.0)*state[P_T][i];
    state[P_et][i]  = state[P_e][i]+0.5*SQR( state[P_ux][i] );
    state[P_ht][i] = state[P_et][i] + state[P_p][i]/state[P_rho][i];
    state[P_h][i]  = state[P_e][i] + state[P_p][i]/state[P_rho][i];
  }
}


void calc_allstates_from_rhopu(double state[NOPRIMS][NOIFX])
{
  int i;
  for (i=0; i< NOIFX; i++)
  {
    state[P_T][i]  = state[P_p][i] / (thermo_R * state[P_rho][i]);
    state[P_e][i]  = thermo_R/(thermo_Gamma-1.0)*state[P_T][i];
    state[P_et][i] = state[P_e][i]+0.5*SQR( state[P_ux][i] );
    state[P_ht][i] = state[P_et][i] + state[P_p][i]/state[P_rho][i];
    state[P_h][i]  = state[P_e][i] + state[P_p][i]/state[P_rho][i];
  }
}


void calc_allP_from_rhopu(double P[NOPRIMS][NOCELLSX])
{
  int i;
  for (i=0; i< NOCELLSX; i++)
  {
    P[P_T][i]  = P[P_p][i] / (thermo_R * P[P_rho][i]);
    P[P_e][i]  = thermo_R/(thermo_Gamma-1.0)*P[P_T][i];
    P[P_et][i] = P[P_e][i]+0.5*SQR( P[P_ux][i] );
    P[P_ht][i] = P[P_et][i] + P[P_p][i]/P[P_rho][i];
    P[P_h][i]  = P[P_e][i] + P[P_p][i]/P[P_rho][i];
  }
}


void calc_allP_from_uetrho(double P[NOPRIMS][NOCELLSX])
{
  int i;
  for (i=0; i < NOCELLSX; i++)
  {
    P[P_e][i]  = P[P_et][i] - 0.5*SQR( P[P_ux][i] );
    P[P_p][i]  = (thermo_Gamma-1.0) * P[P_rho][i] * P[P_e][i];
    P[P_T][i]  = P[P_p][i] / (thermo_R * P[P_rho][i]);
    P[P_ht][i] = P[P_et][i] + P[P_p][i]/P[P_rho][i];
    P[P_h][i]  = P[P_e][i] + P[P_p][i]/P[P_rho][i];
  }//of for
}//of cons_prims()


void calc_allstates_from_rhouht(double states[NOPRIMS][NOIFX])
{
  int i;
  for (i=0; i < NOIFX; i++)
  {
   states[P_h][i]  = states[P_ht][i] - 0.5*SQR( states[P_ux][i] );
   states[P_p][i]  = states[P_h][i] * states[P_rho][i]/(1.0+1.0/(thermo_Gamma-1));
   states[P_e][i]  = states[P_h][i] - states[P_p][i]/states[P_rho][i];
   states[P_et][i] = states[P_e][i] + 0.5*SQR( states[P_ux][i] );
   states[P_T][i]  = states[P_p][i] / (thermo_R * states[P_rho][i]);
  }//of for
}//of cons_prims()



void calc_aM_fromT(double P[NOPRIMS][NOCELLSX], double V[NOADDS][NOCELLSX])
{
  int i;
  for (i=0; i<NOCELLSX; i++)
  {
    V[V_a][i]=get_a(P[P_T][i]);
    V[V_M][i]=P[P_ux][i]/V[V_a][i];
  }//of for
}//of get_a_M()


double get_a(double T)
{
  double a;

  a = sqrt(thermo_Gamma*thermo_R*T);

  return a;
}



