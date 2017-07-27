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

#ifndef THERMODYNAMICS_H
#define THERMODYNAMICS_H

#include "varutils.h"


void calc_allP_from_rhopu(double P[NOPRIMS][NOCELLSX]);

void calc_allstates_from_Tpu(double state[NOPRIMS][NOIFX]);

void calc_allstates_from_rhopu(double state[NOPRIMS][NOIFX]);

void calc_allP_from_Tpu(double P[NOPRIMS][NOCELLSX]);

void calc_allP_from_uetrho(double P[][NOCELLSX]);

void calc_aM_fromT(double P[NOPRIMS][NOCELLSX], double V[NOADDS][NOCELLSX]);

void calc_allstates_from_rhouht(double states[NOPRIMS][NOIFX]);

double get_a(double T);

#endif
