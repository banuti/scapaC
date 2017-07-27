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

#ifndef INICONDITIONS_H
#define INICONDITIONS_H

#include "varutils.h"


void initialconditions(double P[NOEQNS][NOCELLSX], double U[NOEQNS][NOCELLSX]);

void init_shocktube_sod(double P[NOEQNS][NOCELLSX]);

void init_shocktube_Tpu(double P[NOEQNS][NOCELLSX]);

void init_shocktube_rhopu(double P[NOEQNS][NOCELLSX]);

void init_nozzle(double A[NOCELLSX], double dx[NOCELLSX]);

void init_nozzle_Tpu(double P[NOEQNS][NOCELLSX]);

void init_wave_Tpu(double P[NOEQNS][NOCELLSX]);

void init_ambient_Tpu(double P[NOEQNS][NOCELLSX]);

#endif
