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

#ifndef VARUTILS_H
#define VARUTILS_H

#define PI 3.1415926


#define ILENGTH    100
#define IFIRST     1
#define ILAST      (ILENGTH)
#define NOCELLSX   (ILENGTH+2)
#define NOIFX      (NOCELLSX-1)
#define NOEQNS     3
#define NOPRIMS    8
#define NOADDS     2
#define ORDERSPACE 1
#define LIMITER    1




///conservative variables
enum U
{
  U_rho,
  U_rhoux,
  U_rhoet,
};

///convective flux
enum F
{
  F_rhoux   = U_rho,
  F_rhouxux = U_rhoux,
  F_rhoht   = U_rhoet,
};

///primitive variables
enum P
{
  P_rho,
  P_ux,
  P_p,
  P_ht,
  P_e,
  P_et,
  P_T,
  P_h,
};

///additional variables
enum V
{
  V_M, //Mach number
  V_a, //speed of sound
};




void reset_residual(double residual[NOEQNS][NOCELLSX]);



void reset_flux(double flux[NOEQNS][NOIFX]);



void initvariables(double        P[NOPRIMS][NOCELLSX],
                   double        U[NOEQNS][NOCELLSX],
                   double residual[NOEQNS][NOCELLSX],
                   double     flux[NOEQNS][NOIFX],
                   double  l_state[NOPRIMS][NOIFX],
                   double  r_state[NOPRIMS][NOIFX],
                   double       dx[NOCELLSX],
                   double        A[NOCELLSX]);

void fixnegativeprims(double P[NOPRIMS][NOCELLSX]);

void cons_prims(double P[NOPRIMS][NOCELLSX],
                double U[NOEQNS][NOCELLSX]);


void prims_cons(double P[NOPRIMS][NOCELLSX],
                double U[NOEQNS][NOCELLSX]);

double get_timestep(double CFL,
                    double   P[NOPRIMS][NOCELLSX],
                    double  dx[NOCELLSX]);


#endif

