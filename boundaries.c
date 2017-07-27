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

#include <stdio.h>
#include <math.h>

#include "boundaries.h"
#include "thermodynamics.h"
#include "varutils.h"

#define SQR(a)        ( (a)*(a) )


extern const double thermo_Gamma;
extern const double thermo_R;



void update_primcon_bc(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX], double t)
{

  //bc_left_dirichlet(P, U);
  //bc_periodic(P, U);
  //bc_farfield(P, U);
  //bc_nozzle(P, U);
  bc_right_wall(P, U);
  bc_left_wall(P, U);
  //bc_left_piston(P, t);
}


void update_states_bc(double l_state[NOPRIMS][NOIFX],
                      double r_state[NOPRIMS][NOIFX],
                      double t,
                      double P[NOPRIMS][NOCELLSX])
{


  //bc_state_periodic(l_state, r_state);
  //bc_state_farfield(l_state, r_state);
  bc_state_right_wall(l_state, r_state);
  bc_state_left_wall(l_state, r_state);
  //bc_state_left_piston(l_state, r_state, t, P);
}




void update_flux_bc(double P[NOPRIMS][NOCELLSX],
                    double U[NOEQNS][NOCELLSX],
                    double flux[NOEQNS][NOIFX])
{


  bc_flux_right_wall(flux, P);

  bc_flux_left_wall(flux, P);

}




void bc_left_dirichlet(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX])
{

  int    eqn;
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
    
    
    P[P_T][0]     = P[P_T][1]     = T1;
    P[P_p][0]     = P[P_p][1]     = p1;
    P[P_ux][0]    = P[P_ux][1]    = 290.;
    
      
    calc_allP_from_Tpu(P);
  
}


void bc_periodic(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX])
{

  int eqn;
  for (eqn = 0; eqn < NOPRIMS; eqn++)
  {
     P[eqn][0]          = P[eqn][ILAST];
     P[eqn][NOCELLSX-1] = P[eqn][IFIRST];  //NOCELLS-1 is last cell in an array of length NOCELLS
  }
  for (eqn = 0; eqn < NOEQNS; eqn++)
  {
     U[eqn][0]          = U[eqn][ILAST];
     U[eqn][NOCELLSX-1] = U[eqn][IFIRST];
  }

}


void bc_state_periodic(double l_state[NOPRIMS][NOIFX], double r_state[NOPRIMS][NOIFX])
{

  int eqn;
  for (eqn = 0; eqn < NOPRIMS; eqn++)
  {
     l_state[eqn][0]        = l_state[eqn][NOIFX-1];
     r_state[eqn][NOIFX-1]  = r_state[eqn][0];
//NOCELLS-1 is last cell in an array of length NOCELLS
  }

}





void bc_right_wall(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX])
{

 int eqn;
   for (eqn = 0; eqn < NOEQNS; eqn++)
   {
      U[eqn][NOCELLSX-1] = U[eqn][ILAST];
   }//of for


  for (eqn = 0; eqn < NOPRIMS; eqn++)
  {

     P[eqn][NOCELLSX-1] = P[eqn][ILAST];
  }//of for


   U[U_rhoux][NOCELLSX-1] = -U[U_rhoux][ILAST];
  P[P_ux][NOCELLSX-1]    = -P[P_ux][ILAST];


}



void bc_state_right_wall(double l_states[NOPRIMS][NOIFX], double r_states[NOEQNS][NOIFX])
{

 int p;


  for (p = 0; p < NOPRIMS; p++)
  {
     r_states[p][NOIFX-1] = l_states[p][NOIFX-1];

  }//of for

  r_states[P_ux][NOIFX-1] = -l_states[P_ux][NOIFX-1];

}





void bc_flux_left_wall(double flux[NOPRIMS][NOIFX], double P[NOEQNS][NOCELLSX])
{
     flux[F_rhoux][0] = 0.0;
     flux[F_rhouxux][0] = P[P_p][1]-P[P_ux][1]*P[P_rho][1]*sqrt(thermo_Gamma*thermo_R*P[P_T][1]);
     flux[F_rhoht][0] = 0.0;
}




void bc_flux_right_wall(double flux[NOPRIMS][NOIFX], double P[NOEQNS][NOCELLSX])
{
     flux[F_rhoux][NOIFX-1] = 0.0;
     flux[F_rhouxux][NOIFX-1] = P[P_p][1]+P[P_ux][1]*P[P_rho][1]*sqrt(thermo_Gamma*thermo_R*P[P_T][1]);
     flux[F_rhoht][NOIFX-1] = 0.0;
}



void bc_left_wall(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX])
{

 int eqn;
   for (eqn = 0; eqn < NOEQNS; eqn++)
   {
      U[eqn][0] = U[eqn][1];
   }//of for


  for (eqn = 0; eqn < NOPRIMS; eqn++)
  {

     P[eqn][0] = P[eqn][1];
  }//of for


   U[U_rhoux][0] = -U[U_rhoux][1];
  P[P_ux][0]    = -P[P_ux][1];


}


void bc_state_left_wall(double l_states[NOPRIMS][NOIFX], double r_states[NOEQNS][NOIFX])
{

 int p;

  for (p = 0; p < NOPRIMS; p++)
  {

     l_states[p][0] = r_states[p][0];
  }//of for


//   U[U_rhoux][0] = -U[U_rhoux][1];
  l_states[P_ux][0]    = -r_states[P_ux][0];

}








void bc_left_piston(double P[NOPRIMS][NOCELLSX], double t)
{

  double mode=1.0;

//  double k = mode*PI*sqrt(thermo_Gamma*thermo_R*P[P_T][0])/(double)ILENGTH;
  double k = mode*PI*sqrt(thermo_Gamma*thermo_R*0.048780487)/(double)ILENGTH;
//double k = mode*PI*sqrt(thermo_Gamma*thermo_R*P[P_T][NOCELLSX/2])/(double)ILENGTH;

    P[P_ux][0]  = (1.0*sin(k*t))/100.0;

//printf("Piston velocity: u=%f\n", P[P_ux][0]);
    P[P_T][0]   = P[P_T][1];

    P[P_p][0]   = P[P_p][1];



//   U[U_rhoux][0] = -U[U_rhoux][1];

    P[P_rho][0] = P[P_p][0]/(thermo_R * P[P_T][0]);
    P[P_e][0]   = thermo_R/(thermo_Gamma-1.0)*P[P_T][0];
    P[P_et][0]  = P[P_e][0]+0.5*SQR( P[P_ux][0] );
    P[P_ht][0]  = P[P_et][0] + P[P_p][0]/P[P_rho][0];
    P[P_h][0]   = P[P_e][0] + P[P_p][0]/P[P_rho][0];



}


void bc_state_left_piston(double l_states[NOPRIMS][NOIFX],
                          double r_states[NOEQNS][NOIFX],
                          double t, 
                          double P[NOPRIMS][NOCELLSX])
{

 int p;

  for (p = 0; p < NOPRIMS; p++)
  {
     l_states[p][0] = P[p][0];
  }//of for



}







void bc_nozzle(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX])
{
//   P[eqn][0]        = P[eqn][ILAST];
//   P[eqn][NOCELLSX] = P[eqn][IFIRST];

    P[P_T][0]     = 300.0;
    P[P_p][0]     = 1000.0;
    P[P_ux][0]    = 0.0;


/*
    P[P_T][NOCELLSX-1]     = 300.0;
    P[P_p][NOCELLSX-1]     = 10.0;
    P[P_ux][NOCELLSX-1]    = P[P_ux][ILAST];*/


  int eqn;
  for (eqn = 0; eqn < NOEQNS; eqn++)
  {

     U[eqn][NOCELLSX-1] = U[eqn][ILAST];
  }//of for


  for (eqn = 0; eqn < NOPRIMS; eqn++)
  {

     P[eqn][NOCELLSX-1] = P[eqn][ILAST];
  }//of for


}//of bc_periodic

