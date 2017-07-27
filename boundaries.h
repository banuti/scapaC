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

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "varutils.h"



void update_primcon_bc(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX], double t);

void update_flux_bc(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX], double flux[NOEQNS][NOIFX]);

void update_states_bc(double l_state[NOPRIMS][NOIFX],
                      double r_state[NOPRIMS][NOIFX],
                      double t,
                      double P[NOPRIMS][NOCELLSX]);

void bc_periodic(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX]);

void bc_state_periodic(double l_state[NOPRIMS][NOIFX], double r_state[NOPRIMS][NOIFX]);

void bc_farfield(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX]);

void bc_state_farfield(double l_state[NOPRIMS][NOIFX], double r_state[NOPRIMS][NOIFX]);

void bc_left_wall(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX]);

void bc_state_left_wall(double l_states[NOPRIMS][NOIFX], double r_states[NOEQNS][NOIFX]);

void bc_state_right_wall(double l_state[NOPRIMS][NOIFX], double r_state[NOEQNS][NOIFX]);

void bc_right_wall(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX]);

void bc_nozzle(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX]);

void bc_left_piston(double P[NOPRIMS][NOCELLSX], double t);

void bc_state_left_piston(double l_states[NOPRIMS][NOIFX],
                          double r_states[NOEQNS][NOIFX],
                          double t,
                          double P[NOPRIMS][NOCELLSX]);

void bc_flux_right_wall(double flux[NOPRIMS][NOIFX], double P[NOEQNS][NOCELLSX]);

void bc_flux_left_wall(double flux[NOPRIMS][NOIFX], double P[NOEQNS][NOCELLSX]);

void bc_left_dirichlet(double P[NOPRIMS][NOCELLSX], double U[NOEQNS][NOCELLSX]);

#endif
