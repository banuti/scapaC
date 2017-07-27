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
#include <stdio.h>

#include "reconstruction.h"
#include "varutils.h"
#include "thermodynamics.h"


#define MIN(a,b)      ( (a)<=(b) ? (a) : (b) )
#define SQR(a)        ( (a)*(a) )

extern const double thermo_Gamma;
extern const double thermo_R;

void calculate_states(double l_states[NOPRIMS][NOIFX],
                      double r_states[NOPRIMS][NOIFX],
                      double        P[NOPRIMS][NOCELLSX],
                      double       dx[NOCELLSX])
{
  int chooseorder = ORDERSPACE;
  switch(chooseorder)
  {
  case 1:
    first_states(l_states, r_states, P);
    break;

  case 2:
    // second_minmod_states(l_states, r_states, P, dx);
    second_minmod_if_states(l_states, r_states, P, dx);
    //second_allornone_states(l_states, r_states, P, dx);
    break;
  }

}//of calculate states()


void first_states(double l_states[NOPRIMS][NOIFX],
                  double r_states[NOPRIMS][NOIFX],
                  double        P[NOPRIMS][NOCELLSX])
{
  int i,p;
  for (i=0; i<NOIFX; i++)//interface wise loop
  {
    for (p=0; p<NOPRIMS; p++)
    {
      l_states[p][i] = P[p][i];
      r_states[p][i] = P[p][i+1];
    }
  }//of for
}


void second_minmod_states(double l_states[NOPRIMS][NOIFX],
                          double r_states[NOPRIMS][NOIFX],
                          double        P[NOPRIMS][NOCELLSX],
                          double       dx[NOCELLSX])
{
  int i,p;
 // int counterl, counterr, counter0;
  double l_gradient, r_gradient, gradient, rdx;

///cell wise loop
  for (i=IFIRST; i<=ILAST; i++)
  {
//counterl=counterr=counter0=0;
    rdx=1.0/dx[i];
    for (p=0; p<=2; p++)///for rho, ux, ht
    {
///calculate gradients
     l_gradient = (P[p][i]  - P[p][i-1])*rdx;
     r_gradient = (P[p][i+1]- P[p][i]  )*rdx;

    if (LIMITER == 1)
    {
///choose gradient
      if (l_gradient*r_gradient > 0.0) ///if gradients have same sign
      {
       // gradient = MIN((fabs(l_gradient)),(fabs(r_gradient)));
        if (fabs(l_gradient) < fabs(r_gradient))
        {
        gradient = l_gradient;
  // counterl+=1;
  // if ((counterr>0) || (counterr>0))
  //   printf("Inconsistent 2nd Order extrapolation in cell %i\n", i);
        }else
        {
        gradient = r_gradient;
  // counterr+=1;
  // if ((counterl>0) || (counter0>0))
  //   printf("Inconsistent 2nd Order extrapolation in cell %i\n", i);
        }

      }else
      {
        gradient = 0.0;
  // counter0+=1;
  // if ((counterl>0) || (counterr>0))
  //   printf("Inconsistent 2nd Order extrapolation in cell %i\n", i);
      }
    }else
    {
    gradient = 0.5* (l_gradient+r_gradient);
    }

///extrapolate primitive variables (!offset in index!)
    // gradient =0.0;
     r_states[p][i-1] = P[p][i] - gradient * 0.5 * dx[i];
     l_states[p][i]   = P[p][i] + gradient * 0.5 * dx[i];

    }//of for over primitive variables
  }//of for over cells

}



void second_minmod_if_states(double l_states[NOPRIMS][NOIFX],
                          double r_states[NOPRIMS][NOIFX],
                          double        P[NOPRIMS][NOCELLSX],
                          double       dx[NOCELLSX])
{
  int i,p;
 // int counterl, counterr, counter0;
  double l_gradient, r_gradient, gradient, c_gradient, l_grad_lim, r_grad_lim, rdx;

///interface wise loop
  for (i=1; i<NOIFX-1; i++)
  {
//counterl=counterr=counter0=0;
    rdx=1.0/dx[i];
    for (p=0; p<=2; p++)///for rho, ux, ht
    {
///calculate gradients
     l_gradient = (P[p][i]  - P[p][i-1])*rdx;
     c_gradient = (P[p][i+1]- P[p][i]  )*rdx;
     r_gradient = (P[p][i+2]- P[p][i+1])*rdx;


    if (LIMITER == 1)
    {
///choose left gradient
      if (l_gradient*c_gradient > 0.0) ///if gradients have same sign
      {
        if (fabs(l_gradient) < fabs(c_gradient))
         l_grad_lim = l_gradient;
        else
         l_grad_lim = c_gradient;

      }else
        l_grad_lim = 0.0;
///chose right gradient
      if (c_gradient*r_gradient > 0.0) ///if gradients have same sign
      {
        if (fabs(c_gradient) < fabs(r_gradient))
         r_grad_lim = c_gradient;
        else
         r_grad_lim = r_gradient;

      }else
        r_grad_lim = 0.0;


    }else
      r_grad_lim = l_grad_lim = 0.5* (l_gradient+r_gradient);


///extrapolate primitive variables (!offset in index!)
    // gradient =0.0;
     l_states[p][i]  = P[p][i]   + l_grad_lim * 0.5 * dx[i];
     r_states[p][i]  = P[p][i+1] - r_grad_lim * 0.5 * dx[i];


    }//of for over primitive variables
  }//of for over cells


   for (p=0; p<=2; p++)///for rho, ux, ht
    {

     r_states[p][0]       = l_states[p][1]; //P[p][i] - gradient * 0.5 * dx[i];
     l_states[p][NOIFX-1] = r_states[p][NOIFX-2]; //P[p][i] + gradient * 0.5 * dx[i];

    }


}




///macht in dieser Form keinen Sinn, da rho, u und ht sich gerade 
///NICHT gleichsinnig entwickeln - anpassen!
void second_allornone_states(double l_states[NOPRIMS][NOIFX],
                          double r_states[NOPRIMS][NOIFX],
                          double        P[NOPRIMS][NOCELLSX],
                          double       dx[NOCELLSX])
{
  int i,p;
  int counterl, counterr, counter0, allornone;
  double l_gradient[NOPRIMS], r_gradient[NOPRIMS], gradient[NOPRIMS], rdx;

  for (i=IFIRST; i<=ILAST; i++)///cell wise loop
  {
    counterl=counterr=counter0=allornone=0;
    rdx=1.0/dx[i];
    for (p=0; p<=2; p++)///for rho, ux, ht
    {
///calculate gradients
     l_gradient[p] = (P[p][i]  - P[p][i-1])*rdx;
     r_gradient[p] = (P[p][i+1]- P[p][i]  )*rdx;
     }

    for (p=0; p<=2; p++)///for rho, ux, ht
    {
///choose gradient
    if (l_gradient[p]*r_gradient[p] > 0.0) ///if gradients have same sign
    {
     // gradient = MIN((fabs(l_gradient)),(fabs(r_gradient)));
      if (fabs(l_gradient[p]) < fabs(r_gradient[p]))
      {
      gradient[p] = l_gradient[p];
      counterl+=1;
      if ((counterr>0) || (counterr>0))
        allornone+=1;
//   printf("Inconsistent 2nd Order extrapolation in cell %i\n", i);
      }else
      {
      gradient[p] = r_gradient[p];
      counterr+=1;
      if ((counterl>0) || (counter0>0))
        allornone+=1;
//   printf("Inconsistent 2nd Order extrapolation in cell %i\n", i);
      }

    }else
    {
      gradient[p] = 0.0;
      counter0+=1;
      if ((counterl>0) || (counterr>0))
        allornone+=1;
//   printf("Inconsistent 2nd Order extrapolation in cell %i\n", i);
    }

    if (allornone !=0)
    {
      for (p=0; p<=2; p++)///for rho, ux, ht
      {
         gradient[p] = 0.0;
      }
    }

    for (p=0; p<=2; p++)///for rho, ux, ht
    {
///extrapolate primitive variables (!offset in index!)
    // gradient =0.0;
     //gradient = 0.5* (l_gradient+r_gradient);
     r_states[p][i-1] = P[p][i] - gradient[p] * 0.5 * dx[i];
     l_states[p][i]   = P[p][i] + gradient[p] * 0.5 * dx[i];

    }//of for over primitive variables
  }//of for over cells



  for (p=0; p<NOPRIMS; p++)
  {
    l_states[p][0]       = P[p][0];
    r_states[p][NOIFX-1] = P[p][NOCELLSX-1];
  }

  calc_allstates_from_rhouht(l_states);
  calc_allstates_from_rhouht(r_states);


  }//of cell loop

}//of function


