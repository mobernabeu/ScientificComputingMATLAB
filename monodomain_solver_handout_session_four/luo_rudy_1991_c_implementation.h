// Introduction to Scientific and High Performance Computing with MATLAB
//
// Copyright (c) 2012-2014, Miguel O. Bernabeu, All rights reserved.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library.

//==============================================================================
// CellML file:   C:\Users\Miguel Bernabeu\Desktop\luo_rudy_1991.cellml
// CellML model:  luo_rudy_1991
// Date and time: 4/13/2012 at 11:30:15 at AM
//------------------------------------------------------------------------------
// Conversion from CellML 1.0 to C (header) was done using COR (0.9.31.1409)
//    Copyright 2002-2012 Dr Alan Garny
//    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
//------------------------------------------------------------------------------
// http://www.cellml.org/
//==============================================================================

#ifndef __LUO_RUDY_1991_H__
#define __LUO_RUDY_1991_H__

//------------------------------------------------------------------------------
// State variables
//------------------------------------------------------------------------------

#define _NB_OF_STATE_VARIABLES_ 8

//extern double Y[_NB_OF_STATE_VARIABLES_];
//extern double dY[_NB_OF_STATE_VARIABLES_];
// 0: h (dimensionless) (in fast_sodium_current_h_gate)
// 1: j (dimensionless) (in fast_sodium_current_j_gate)
// 2: m (dimensionless) (in fast_sodium_current_m_gate)
// 3: Cai (millimolar) (in intracellular_calcium_concentration)
// 4: V (millivolt) (in membrane)
// 5: d (dimensionless) (in slow_inward_current_d_gate)
// 6: f (dimensionless) (in slow_inward_current_f_gate)
// 7: X (dimensionless) (in time_dependent_potassium_current_X_gate)

extern char YNames[_NB_OF_STATE_VARIABLES_][4];
extern char YUnits[_NB_OF_STATE_VARIABLES_][14];
extern char YComponents[_NB_OF_STATE_VARIABLES_][40];

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

extern double E_b;   // millivolt (in background_current)
extern double g_b;   // milliS_per_cm2 (in background_current)
extern double g_Na;   // milliS_per_cm2 (in fast_sodium_current)
extern double Ki;   // millimolar (in ionic_concentrations)
extern double Ko;   // millimolar (in ionic_concentrations)
extern double Nai;   // millimolar (in ionic_concentrations)
extern double Nao;   // millimolar (in ionic_concentrations)
extern double C;   // microF_per_cm2 (in membrane)
extern double F;   // coulomb_per_mole (in membrane)
extern double R;   // joule_per_kilomole_kelvin (in membrane)
extern double T;   // kelvin (in membrane)
extern double stim_amplitude;   // microA_per_cm2 (in membrane)
extern double stim_duration;   // millisecond (in membrane)
extern double stim_end;   // millisecond (in membrane)
extern double stim_period;   // millisecond (in membrane)
extern double stim_start;   // millisecond (in membrane)
extern double g_Kp;   // milliS_per_cm2 (in plateau_potassium_current)
extern double PR_NaK;   // dimensionless (in time_dependent_potassium_current)

//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------

extern double i_b;   // microA_per_cm2 (in background_current)
extern double alpha_h;   // per_millisecond (in fast_sodium_current_h_gate)
extern double beta_h;   // per_millisecond (in fast_sodium_current_h_gate)
extern double alpha_j;   // per_millisecond (in fast_sodium_current_j_gate)
extern double beta_j;   // per_millisecond (in fast_sodium_current_j_gate)
extern double alpha_m;   // per_millisecond (in fast_sodium_current_m_gate)
extern double beta_m;   // per_millisecond (in fast_sodium_current_m_gate)
extern double E_Na;   // millivolt (in fast_sodium_current)
extern double i_Na;   // microA_per_cm2 (in fast_sodium_current)
extern double I_stim;   // microA_per_cm2 (in membrane)
extern double E_Kp;   // millivolt (in plateau_potassium_current)
extern double Kp;   // dimensionless (in plateau_potassium_current)
extern double i_Kp;   // microA_per_cm2 (in plateau_potassium_current)
extern double alpha_d;   // per_millisecond (in slow_inward_current_d_gate)
extern double beta_d;   // per_millisecond (in slow_inward_current_d_gate)
extern double alpha_f;   // per_millisecond (in slow_inward_current_f_gate)
extern double beta_f;   // per_millisecond (in slow_inward_current_f_gate)
extern double E_si;   // millivolt (in slow_inward_current)
extern double i_si;   // microA_per_cm2 (in slow_inward_current)
extern double alpha_X;   // per_millisecond (in time_dependent_potassium_current_X_gate)
extern double beta_X;   // per_millisecond (in time_dependent_potassium_current_X_gate)
extern double Xi;   // dimensionless (in time_dependent_potassium_current_Xi_gate)
extern double E_K;   // millivolt (in time_dependent_potassium_current)
extern double g_K;   // milliS_per_cm2 (in time_dependent_potassium_current)
extern double i_K;   // microA_per_cm2 (in time_dependent_potassium_current)
extern double K1_infinity;   // dimensionless (in time_independent_potassium_current_K1_gate)
extern double alpha_K1;   // per_millisecond (in time_independent_potassium_current_K1_gate)
extern double beta_K1;   // per_millisecond (in time_independent_potassium_current_K1_gate)
extern double E_K1;   // millivolt (in time_independent_potassium_current)
extern double g_K1;   // milliS_per_cm2 (in time_independent_potassium_current)
extern double i_K1;   // microA_per_cm2 (in time_independent_potassium_current)

extern int cellModelInitialised;

//------------------------------------------------------------------------------
// Procedures
//------------------------------------------------------------------------------

extern void init();
extern void ComputeTimeDerivative(double time, double *previousSolution, double *computedTimeDerivative, double stimulusCurrent);
extern double ComputeIonicCurrent(double *previousSolution);

//------------------------------------------------------------------------------

#endif

//==============================================================================
// End of file
//==============================================================================
