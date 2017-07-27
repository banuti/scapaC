/***************************************************************************
 *   Copyright (C) 2008 by Daniel Banuti   *
 *   daniel@banuti.de   *
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "varutils.h"
#include "boundaries.h"
#include "ioutils.h"
#include "reconstruction.h"
#include "thermodynamics.h"
#include "fluxes.h"
#include "iniconditions.h"




////read param


int main(int argc, char *argv[])
{

//init consts

// double     initdt      = 0.0000001;
double     tmax        = 0.15;
int        nmax        = 600;
int        outinterval = 50;
double     CFL         = 0.8;  // original Sod: 0.9, Liou1996: 0.8
double     dt          = 0.0000001;
double     t           = 0.0;

double     dx[NOCELLSX];
double      A[NOCELLSX];

////double gradients[nocells];

//typedef cell
double            P[NOPRIMS][NOCELLSX];
double            U[NOEQNS][NOCELLSX];
double            V[NOADDS][NOCELLSX];
double            residual[NOEQNS][NOCELLSX];
double            flux[NOEQNS][NOIFX];
double            l_state[NOPRIMS][NOIFX];
double            r_state[NOPRIMS][NOIFX];





 printf("\n Hallo Welt!\n\n");


///sets fields, variables to zero
  initvariables(P, U, residual, flux, l_state, r_state, dx, A);


// outputtecplotheader();
  outputgnuplotheader();
  outhistoryheader();

///init vecUfield
  initialconditions(P, U);
  

  //init_nozzle(A, dx);

  prims_cons(P, U);


///output values
  outputscreen(-1, P, A);
  outputfile(-1, P, A, t);

///TIMELOOP
  int n;
  for (n=0; n<=nmax; n++)
  {
//outputfile(-85, P, A);
///timestepping

printf("\n Hallo Welt before!\n\n");

  dt = get_timestep(CFL, P, dx);
  // printf('dt %f',dt);
  t=dt*(float)n;

///fill ghost cells acc to bc
    update_primcon_bc(P, U, t);
//outputfile(-80, P, A);
//outputscreen(-10, P, A);

///reconstruction of states
    calculate_states(l_state, r_state, P, dx);
//outputfile(-70, P, A);
// outputscreen(-20, P, A);
//outputscreen(-20, l_state, A);
//outputscreen(-20, r_state, A);
///change states according to bc
   update_states_bc(l_state,r_state, t, P);
// outputscreen(-30, P, A);

///calculate all state variables from rho, u, ht reconstruction
//   calc_allstates_from_rhouht(l_state);
//   calc_allstates_from_rhouht(r_state);

calc_allstates_from_rhopu(l_state);
calc_allstates_from_rhopu(r_state);

// outputscreen(-40, P, A);
//outputscreen(-40, l_state, A);
//outputscreen(-40, r_state, A);
///flux determination
    compute_fluxes(l_state, r_state, flux, A);
//outputfile(-60, P, A);
//  outputscreen(-50, P, A);

///change fluxes according to bc
//    update_flux_bc(P, U, flux);
//outputscreen(-60, P, A);
///update residuals
    update_residuals(flux, residual, U, dt, dx, A);
    reset_flux(flux);
    reset_residual(residual);
//outputfile(-50, P, A);
//outputscreen(-70, P, A);
  ///cons-prims
    cons_prims(P, U);
    
//outputscreen(-75, P, A);
    calc_allP_from_uetrho(P);
    
    
    fixnegativeprims(P);
//    calc_aM_fromT(P, V);
//  outputscreen(-80, P, A);
  ///output
    if (fmod(n,outinterval) == 0)
    {
      printf("t=%f, dt=%f\n", t, dt);
      // outputscreen(n, P, A);
      outputfile(n, P, A, t);
    }

  }///of timeloop

  //output values
  printf("t=%f, dt=%f\n", t, dt);
  outputscreen(n, P, A);
  outputfile(n, P, A, t);

  return EXIT_SUCCESS;
}
