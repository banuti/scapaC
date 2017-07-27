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

#include "ioutils.h"
#include "varutils.h"
#include "thermodynamics.h"

extern const double thermo_Gamma;
extern const double thermo_R;


void outputscreen(int n, double P[NOPRIMS][NOCELLSX], double A[NOCELLSX])
{
  printf("=========== Timestep %i =============\n",n);
  printf("cell\t u\t\t rho\t\t p\t\t T\t\t A\t\t M\t\t mdot\n");
  int i;
  double M, mdot;
  for (i=0; i<NOCELLSX; i++)
  {
    M=P[P_ux][i]/sqrt(thermo_Gamma*thermo_R*P[P_T][i]);
    mdot=P[P_rho][i]*P[P_ux][i]*A[i];
    printf("%i\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",i,P[P_ux][i],P[P_rho][i],P[P_p][i],P[P_T][i], A[i], M, mdot);
  }
}



void outputfile(int n, double P[NOPRIMS][NOCELLSX], double A[NOCELLSX], double t)
{

//  outputfield(n, P, A);
  outputgnufield(n, P, A);
  outputhistory(n, P, t);

}



void outputtecplotheader(void)
{
  FILE *file_ptr;
  file_ptr = fopen("outfield.dat","w");
  fclose(file_ptr);
  file_ptr = fopen("outfield.dat","a");
  fprintf(file_ptr,"VARIABLES = \"cell\" \"u\" \"rho\" \"p\" \"T\" \"A\" \n");
  fclose(file_ptr);
}

void outputfield(int n, double P[NOPRIMS][NOCELLSX], double A[NOCELLSX])
{

  printf("Write output\n");

  FILE *file_ptr;
  file_ptr = fopen("outfield.dat","a");

  fprintf(file_ptr,"ZONE T=\"Timestep %i\" \n",n);

  int i;
  for (i=0; i<NOCELLSX; i++)
  {
    fprintf(file_ptr,"%i\t %f\t %f\t %f\t %f\t %f\n",i,P[P_ux][i],P[P_rho][i],P[P_p][i],P[P_T][i],A[i]);
  }

  fclose(file_ptr);
}


void outputgnuplotheader(void)
{
    FILE *file_ptr;
    file_ptr = fopen("outfield.dat","w");
    fclose(file_ptr);
    file_ptr = fopen("outfield.dat","a");
    fprintf(file_ptr,"# \"timestep\" \"cell\" \"u\" \"rho\" \"p\" \"T\" \"A\" \n");
    fclose(file_ptr);
}

void outputgnufield(int n, double P[NOPRIMS][NOCELLSX], double A[NOCELLSX])
{
    
    printf("Write output field gnuplot\n");
    
    FILE *file_ptr;
    file_ptr = fopen("outfield.dat","a");
        
    int i;
    for (i=0; i<NOCELLSX; i++)
    {
        fprintf(file_ptr,"%i\t %i\t %f\t %f\t %f\t %f\t %f\n",n,i,P[P_ux][i],P[P_rho][i],P[P_p][i],P[P_T][i],A[i]);
    }
    fprintf(file_ptr,"\n\n");
    fclose(file_ptr);
}



void outhistoryheader(void)
{
  FILE *histfile_ptr;
  histfile_ptr = fopen("outhistory.dat","w");
  fclose(histfile_ptr);
  histfile_ptr = fopen("outhistory.dat","a");
  fprintf(histfile_ptr,"VARIABLES = \"t\" \"p_mean\" \"p_end\" \"p_mid\" \"T_mean\" \"m_ges\"\n");
  fclose(histfile_ptr);
}





void outputhistory(int n, double P[NOPRIMS][NOCELLSX], double t)
{

  double p_mean, p_mid, p_end, T_mean, m_ges;

  p_mean = p_mid = p_end = T_mean = m_ges = 0.0;

  printf("Write output history\n");

  FILE *file_ptr;
  file_ptr = fopen("outhistory.dat","a");

  int i;
  for (i=IFIRST; i<=ILAST; i++)
  {
    p_mean += P[P_p][i];
    T_mean += P[P_T][i];
    m_ges  += P[P_rho][i]; //for equidistant dx and constant A
  }

  p_mean /= ILENGTH;
  T_mean /= ILENGTH;


  fprintf(file_ptr,"%f\t %f\t %f\t %f\t %f\t %f\n",t,p_mean,P[P_p][NOCELLSX-1],P[P_p][NOCELLSX/2],P[P_T][i], m_ges);

  fclose(file_ptr);
}


