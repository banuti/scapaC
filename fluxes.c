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

#include "fluxes.h"

#include "thermodynamics.h"


#define SIGN(a)       ( ((a) >= 0.0) ? 1.0 : -1.0)
#define GTZERO(a,b,c) ( ((a) >= 0.0) ? (b) : (c) )
#define SQR(a)        ( (a)*(a) ) 


extern const double thermo_Gamma;
extern const double thermo_R;




void compute_fluxes(double l_states[NOPRIMS][NOIFX],
                    double r_states[NOPRIMS][NOIFX],
                    double flux[NOEQNS][NOIFX],
                    double A[NOCELLSX])
{

  inv_fluxes(l_states, r_states, flux, A);
  visc_fluxes(l_states, r_states, flux);

}



void inv_fluxes(double l_states[NOPRIMS][NOIFX],
                double r_states[NOPRIMS][NOIFX],
                double flux[NOEQNS][NOIFX],
                double A[NOCELLSX])
{
  double ausmp_flux[NOPRIMS];
  double one_r_state[NOPRIMS];
  double one_l_state[NOPRIMS];
  double areaif;

  int i,p;
  for (i=0; i<NOIFX; i++)
  {
    for (p=0; p<NOPRIMS; p++)
    {
      one_r_state[p] = r_states[p][i];
      one_l_state[p] = l_states[p][i];
    }/*
printf("Interface %i \n",i);*/
    areaif = 0.5 * (A[i]+A[i+1]);
    //areaif = A[i];
    ausmp(one_l_state, one_r_state, ausmp_flux, areaif);

    for (p=0; p<NOEQNS; p++)
    {
      //flux[p][i] += areaif * ausmp_flux[p];
      flux[p][i] += ausmp_flux[p];
    }

  }
}






void ausmp(double one_l_state[NOPRIMS],
           double one_r_state[NOPRIMS],
           double ausmp_flux[NOEQNS],
           double areaif)
{
  double var_l[NOEQNS];
  double var_r[NOEQNS];

  double rho_l, vx_l, p_l, a_l, ht_l,
         rho_r, vx_r, p_r, a_r, ht_r;

  double rpc, rrhoc;

  double at_sqr_l, astar_sqr_l, astar_l, atild_l, M_l, Mp, Pp,
         at_sqr_r, astar_sqr_r, astar_r, atild_r, M_r, Mm, Pm, a_half;


  double m_half, m_half_p, m_half_m, p_half;

///AUSM+ constants
  const double alpha = 3.0/16.0; 
  const double beta  = 1.0/8.0;

// ///Van Leer:
  // const double alpha = 0.0; 
  // const double beta  = 0.0;

  rho_l = one_l_state[P_rho];
  rho_r = one_r_state[P_rho];
  vx_l  = one_l_state[P_ux];
  vx_r  = one_r_state[P_ux];
  p_l   = one_l_state[P_p];
  p_r   = one_r_state[P_p];
  ht_l  = one_l_state[P_ht];
  ht_r  = one_r_state[P_ht];
/*
printf("        L \t\t| R \n");
printf("---------------------------\n");
printf(" rho  %f \t\t| %f\n",rho_l,rho_r);
printf(" ux   %f \t\t| %f\n",vx_l,vx_r);
printf(" p    %f \t| %f\n",p_l,p_r);
printf(" ht   %f \t| %f\n",ht_l,ht_r);
printf("\n");*/


///common speed of sound at interface

  //at_sqr_l    = (thermo_Gamma-1)*one_l_state[P_h];
  //at_sqr_r    = (thermo_Gamma-1)*one_r_state[P_h];



  // astar_sqr_l = 2.0*(thermo_Gamma-1)/(thermo_Gamma+1)*one_l_state[P_h];
  // astar_sqr_r = 2.0*(thermo_Gamma-1)/(thermo_Gamma+1)*one_r_state[P_h];

  astar_sqr_l = 2.0*(thermo_Gamma-1)/(thermo_Gamma+1)*ht_l;
  astar_sqr_r = 2.0*(thermo_Gamma-1)/(thermo_Gamma+1)*ht_r;




  astar_l = sqrt(astar_sqr_l);
  astar_r = sqrt(astar_sqr_r);


  atild_l = astar_sqr_l / fmax((astar_l),(fabs(vx_l)));
  atild_r = astar_sqr_r / fmax((astar_r),(fabs(vx_r)));



/// Different possibilities to define a common interface speed of sound a_half
  // a_half = fmin((atild_l),(atild_r));

  a_half = 0.5*(sqrt(thermo_Gamma*thermo_R*one_l_state[P_T])+sqrt(thermo_Gamma*thermo_R*one_r_state[P_T]));
  // a_half = 0.5*(sqrt(thermo_Gamma*one_l_state[P_p]/one_l_state[P_rho])
               // +sqrt(thermo_Gamma*one_r_state[P_p]/one_r_state[P_rho]));

///Mach number based on common speed of sound
  M_l = vx_l/a_half;
  M_r = vx_r/a_half;


///splitting
  if (fabs(M_l) >= 1.0)
  {
    Mp = 0.5*(M_l+fabs(M_l));
    Pp = 0.5*(1+SIGN(M_l));
  } else
  {
    Mp = 0.25*SQR((M_l+1))+beta*SQR((M_l*M_l-1));
    Pp = 0.25*SQR((M_l+1))*(2.0-M_l)+alpha*M_l*SQR((M_l*M_l-1));
  }

  if (fabs(M_r) >= 1.0)
  {
    Mm = 0.5*(M_r-fabs(M_r));
    Pm = 0.5*(1-SIGN(M_r));
  } else
  {
    Mm = -0.25*SQR((M_r-1))-beta*SQR((M_r*M_r-1));
    Pm = 0.25*SQR((M_r-1))*(2.0+M_r)-alpha*M_l*SQR((M_r*M_r-1));
  }


/// M and P at interface
  m_half   = Mp + Mm;
  m_half_p = 0.5*(m_half+fabs(m_half));
  m_half_m = 0.5*(m_half-fabs(m_half));

  p_half = Pp*p_l + Pm*p_r;

/*
printf("Mm    =%f\n",Mm);
printf("Mp    =%f\n",Mp);
printf("m_half=%f\n",m_half);
printf("Pm    =%f\n",Pm);
printf("Pp    =%f\n",Pp);
printf("p_half=%f\n",p_half);*/


///flux
//   ausmp_flux[F_rhoux]   = areaif * ( a_half * ( m_half_p*rho_l + m_half_m*rho_r)); 
//   ausmp_flux[F_rhouxux] = areaif * ( a_half * ( m_half_p*rho_l*vx_l + m_half_m*rho_r*vx_r) + p_half);
//   ausmp_flux[F_rhoht]   = areaif * ( a_half * ( m_half_p*rho_l*ht_l + m_half_m*rho_r*ht_r)); 

  ausmp_flux[F_rhoux]   = ( a_half * ( m_half_p*rho_l + m_half_m*rho_r)); 
  ausmp_flux[F_rhouxux] = ( a_half * ( m_half_p*rho_l*vx_l + m_half_m*rho_r*vx_r) + p_half);
  ausmp_flux[F_rhoht]   = ( a_half * ( m_half_p*rho_l*ht_l + m_half_m*rho_r*ht_r));

/*
printf("Flux F_rhoux   %f\n",ausmp_flux[F_rhoux]);
printf("Flux F_rhouxux %f\n",ausmp_flux[F_rhouxux]);
printf("Flux F_rhoht   %f\n\n",ausmp_flux[F_rhoht]);*/


} /** flux_ausmp() **/








void visc_fluxes(double l_states[NOPRIMS][NOIFX],
                 double r_states[NOPRIMS][NOIFX],
                 double flux[NOEQNS][NOIFX])
{
  ///empty
}




void update_residuals(double     flux[NOEQNS][NOIFX],
                      double residual[NOEQNS][NOCELLSX],
                      double        U[NOEQNS][NOCELLSX],
                      double dt,
                      double dx[NOCELLSX],
                      double A[NOCELLSX])
{

  int i,eqn;

  for (i=IFIRST; i<=ILAST; i++)
  {
//printf("Cell %i\n",i);
    for (eqn=0; eqn<NOEQNS; eqn++)
    {
      residual[eqn][i] += flux[eqn][i-1];
      residual[eqn][i] -= flux[eqn][i];
/*      residual[eqn][i] += 0.5*(A[i-1]+A[i]) * flux[eqn][i-1];
        residual[eqn][i] -= 0.5*(A[i]+A[i+1]) * flux[eqn][i];*/
/*printf("residual %f\n",residual[eqn][i]);*/

    }

  }

  for (i=IFIRST; i<=ILAST; i++)
  {
    for (eqn=0; eqn<NOEQNS; eqn++)
    {
//      U[eqn][i] += dt/( dx[i]*A[i] )*residual[eqn][i];
      U[eqn][i] += dt/dx[i]*residual[eqn][i];
     //U[eqn][i] += A[i]*dt/dx[i]*residual[eqn][i];
    }
  }

}

