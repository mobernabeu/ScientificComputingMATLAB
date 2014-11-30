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
// Conversion from CellML 1.0 to C was done using COR (0.9.31.1409)
//    Copyright 2002-2012 Dr Alan Garny
//    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
//------------------------------------------------------------------------------
// http://www.cellml.org/
//==============================================================================

#include "luo_rudy_1991_c_implementation.h"
#include "mex.h"

//------------------------------------------------------------------------------

#include <math.h>
#include <string.h>

//------------------------------------------------------------------------------
// State variables
//------------------------------------------------------------------------------

//double Y[_NB_OF_STATE_VARIABLES_];
//double dY[_NB_OF_STATE_VARIABLES_];
// 0: h (dimensionless) (in fast_sodium_current_h_gate)
// 1: j (dimensionless) (in fast_sodium_current_j_gate)
// 2: m (dimensionless) (in fast_sodium_current_m_gate)
// 3: Cai (millimolar) (in intracellular_calcium_concentration)
// 4: V (millivolt) (in membrane)
// 5: d (dimensionless) (in slow_inward_current_d_gate)
// 6: f (dimensionless) (in slow_inward_current_f_gate)
// 7: X (dimensionless) (in time_dependent_potassium_current_X_gate)

char YNames[_NB_OF_STATE_VARIABLES_][4];
char YUnits[_NB_OF_STATE_VARIABLES_][14];
char YComponents[_NB_OF_STATE_VARIABLES_][40];

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

double E_b;   // millivolt (in background_current)
double g_b;   // milliS_per_cm2 (in background_current)
double g_Na;   // milliS_per_cm2 (in fast_sodium_current)
double Ki;   // millimolar (in ionic_concentrations)
double Ko;   // millimolar (in ionic_concentrations)
double Nai;   // millimolar (in ionic_concentrations)
double Nao;   // millimolar (in ionic_concentrations)
double C;   // microF_per_cm2 (in membrane)
double F;   // coulomb_per_mole (in membrane)
double R;   // joule_per_kilomole_kelvin (in membrane)
double T;   // kelvin (in membrane)
double stim_amplitude;   // microA_per_cm2 (in membrane)
double stim_duration;   // millisecond (in membrane)
double stim_end;   // millisecond (in membrane)
double stim_period;   // millisecond (in membrane)
double stim_start;   // millisecond (in membrane)
double g_Kp;   // milliS_per_cm2 (in plateau_potassium_current)
double PR_NaK;   // dimensionless (in time_dependent_potassium_current)

//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------

double i_b;   // microA_per_cm2 (in background_current)
double alpha_h;   // per_millisecond (in fast_sodium_current_h_gate)
double beta_h;   // per_millisecond (in fast_sodium_current_h_gate)
double alpha_j;   // per_millisecond (in fast_sodium_current_j_gate)
double beta_j;   // per_millisecond (in fast_sodium_current_j_gate)
double alpha_m;   // per_millisecond (in fast_sodium_current_m_gate)
double beta_m;   // per_millisecond (in fast_sodium_current_m_gate)
double E_Na;   // millivolt (in fast_sodium_current)
double i_Na;   // microA_per_cm2 (in fast_sodium_current)
double I_stim;   // microA_per_cm2 (in membrane)
double E_Kp;   // millivolt (in plateau_potassium_current)
double Kp;   // dimensionless (in plateau_potassium_current)
double i_Kp;   // microA_per_cm2 (in plateau_potassium_current)
double alpha_d;   // per_millisecond (in slow_inward_current_d_gate)
double beta_d;   // per_millisecond (in slow_inward_current_d_gate)
double alpha_f;   // per_millisecond (in slow_inward_current_f_gate)
double beta_f;   // per_millisecond (in slow_inward_current_f_gate)
double E_si;   // millivolt (in slow_inward_current)
double i_si;   // microA_per_cm2 (in slow_inward_current)
double alpha_X;   // per_millisecond (in time_dependent_potassium_current_X_gate)
double beta_X;   // per_millisecond (in time_dependent_potassium_current_X_gate)
double Xi;   // dimensionless (in time_dependent_potassium_current_Xi_gate)
double E_K;   // millivolt (in time_dependent_potassium_current)
double g_K;   // milliS_per_cm2 (in time_dependent_potassium_current)
double i_K;   // microA_per_cm2 (in time_dependent_potassium_current)
double K1_infinity;   // dimensionless (in time_independent_potassium_current_K1_gate)
double alpha_K1;   // per_millisecond (in time_independent_potassium_current_K1_gate)
double beta_K1;   // per_millisecond (in time_independent_potassium_current_K1_gate)
double E_K1;   // millivolt (in time_independent_potassium_current)
double g_K1;   // milliS_per_cm2 (in time_independent_potassium_current)
double i_K1;   // microA_per_cm2 (in time_independent_potassium_current)

int cellModelInitialised = 0;

//------------------------------------------------------------------------------
// Initialisation
//------------------------------------------------------------------------------

void init()
{
   //---------------------------------------------------------------------------
   // State variables
   //---------------------------------------------------------------------------

//    Y[0] = 0.982660523699656;   // h (dimensionless) (in fast_sodium_current_h_gate)
//    Y[1] = 0.989108212766685;   // j (dimensionless) (in fast_sodium_current_j_gate)
//    Y[2] = 0.00171338077730188;   // m (dimensionless) (in fast_sodium_current_m_gate)
//    Y[3] = 0.00017948816388306;   // Cai (millimolar) (in intracellular_calcium_concentration)
//    Y[4] = -84.3801107371;   // V (millivolt) (in membrane)
//    Y[5] = 0.00302126301779861;   // d (dimensionless) (in slow_inward_current_d_gate)
//    Y[6] = 0.999967936476325;   // f (dimensionless) (in slow_inward_current_f_gate)
//    Y[7] = 0.0417603108167287;   // X (dimensionless) (in time_dependent_potassium_current_X_gate)
    
   strcpy(YNames[0], "h");
   strcpy(YNames[1], "j");
   strcpy(YNames[2], "m");
   strcpy(YNames[3], "Cai");
   strcpy(YNames[4], "V");
   strcpy(YNames[5], "d");
   strcpy(YNames[6], "f");
   strcpy(YNames[7], "X");

   strcpy(YUnits[0], "dimensionless");
   strcpy(YUnits[1], "dimensionless");
   strcpy(YUnits[2], "dimensionless");
   strcpy(YUnits[3], "millimolar");
   strcpy(YUnits[4], "millivolt");
   strcpy(YUnits[5], "dimensionless");
   strcpy(YUnits[6], "dimensionless");
   strcpy(YUnits[7], "dimensionless");

   strcpy(YComponents[0], "fast_sodium_current_h_gate");
   strcpy(YComponents[1], "fast_sodium_current_j_gate");
   strcpy(YComponents[2], "fast_sodium_current_m_gate");
   strcpy(YComponents[3], "intracellular_calcium_concentration");
   strcpy(YComponents[4], "membrane");
   strcpy(YComponents[5], "slow_inward_current_d_gate");
   strcpy(YComponents[6], "slow_inward_current_f_gate");
   strcpy(YComponents[7], "time_dependent_potassium_current_X_gate");

   //---------------------------------------------------------------------------
   // Constants
   //---------------------------------------------------------------------------

   E_b = -59.87;   // millivolt (in background_current)
   g_b = 0.03921;   // milliS_per_cm2 (in background_current)
   g_Na = 23.0;   // milliS_per_cm2 (in fast_sodium_current)
   Ki = 145.0;   // millimolar (in ionic_concentrations)
   Ko = 5.4;   // millimolar (in ionic_concentrations)
   Nai = 18.0;   // millimolar (in ionic_concentrations)
   Nao = 140.0;   // millimolar (in ionic_concentrations)
   C = 1.0;   // microF_per_cm2 (in membrane)
   F = 96484.6;   // coulomb_per_mole (in membrane)
   R = 8314.0;   // joule_per_kilomole_kelvin (in membrane)
   T = 310.0;   // kelvin (in membrane)
   stim_amplitude = -25.5;   // microA_per_cm2 (in membrane)
   stim_duration = 2.0;   // millisecond (in membrane)
   stim_end = 9000.0;   // millisecond (in membrane)
   stim_period = 1000.0;   // millisecond (in membrane)
   stim_start = 100.0;   // millisecond (in membrane)
   g_Kp = 0.0183;   // milliS_per_cm2 (in plateau_potassium_current)
   PR_NaK = 0.01833;   // dimensionless (in time_dependent_potassium_current)

   //---------------------------------------------------------------------------
   // Computed variables
   //---------------------------------------------------------------------------

   E_Na = R*T/F*log(Nao/Nai);
   g_K = 0.282*sqrt(Ko/5.4);
   E_K = R*T/F*log((Ko+PR_NaK*Nao)/(Ki+PR_NaK*Nai));
   g_K1 = 0.6047*sqrt(Ko/5.4);
   E_K1 = R*T/F*log(Ko/Ki);
   E_Kp = E_K1;
   
   cellModelInitialised = 1;
}

//------------------------------------------------------------------------------
// Computation
//------------------------------------------------------------------------------

void ComputeTimeDerivative(double time, double *previousSolution, double *computedTimeDerivative, double stimulusCurrent)
{
   // time: time (millisecond)

   if (!cellModelInitialised)
   {
       mexErrMsgIdAndTxt( "MATLAB:luo_rudy:not_initialised",
                "init() must be called before ComputeTimeDerivative");
   }
    
   if (previousSolution==NULL || computedTimeDerivative==NULL)
   {
       mexErrMsgIdAndTxt( "MATLAB:luo_rudy:null_pointer",
                "Either previousSolution or computedTimeDerivative don't point at a valid memory segment");
   }
    
   i_b = g_b*(previousSolution[4]-E_b);
   i_Na = g_Na*pow(previousSolution[2], 3.0)*previousSolution[0]*previousSolution[1]*(previousSolution[4]-E_Na);

   if (previousSolution[4] < -40.0)
      alpha_h = 0.135*exp((80.0+previousSolution[4])/-6.8);
   else
      alpha_h = 0.0;

   if (previousSolution[4] < -40.0)
      beta_h = 3.56*exp(0.079*previousSolution[4])+310000.0*exp(0.35*previousSolution[4]);
   else
      beta_h = 1.0/(0.13*(1.0+exp((previousSolution[4]+10.66)/-11.1)));

   computedTimeDerivative[0] = alpha_h*(1.0-previousSolution[0])-beta_h*previousSolution[0];

   if (previousSolution[4] < -40.0)
      alpha_j = (-127140.0*exp(0.2444*previousSolution[4])-0.00003474*exp(-0.04391*previousSolution[4]))*(previousSolution[4]+37.78)/(1.0+exp(0.311*(previousSolution[4]+79.23)));
   else
      alpha_j = 0.0;

   if (previousSolution[4] < -40.0)
      beta_j = 0.1212*exp(-0.01052*previousSolution[4])/(1.0+exp(-0.1378*(previousSolution[4]+40.14)));
   else
      beta_j = 0.3*exp(-0.0000002535*previousSolution[4])/(1.0+exp(-0.1*(previousSolution[4]+32.0)));

   computedTimeDerivative[1] = alpha_j*(1.0-previousSolution[1])-beta_j*previousSolution[1];
   alpha_m = 0.32*(previousSolution[4]+47.13)/(1.0-exp(-0.1*(previousSolution[4]+47.13)));
   beta_m = 0.08*exp(-previousSolution[4]/11.0);
   computedTimeDerivative[2] = alpha_m*(1.0-previousSolution[2])-beta_m*previousSolution[2];
   E_si = 7.7-13.0287*log(previousSolution[3]/1.0);
   i_si = 0.09*previousSolution[5]*previousSolution[6]*(previousSolution[4]-E_si);
   computedTimeDerivative[3] = -0.0001/1.0*i_si+0.07*(0.0001-previousSolution[3]);

   if ((time >= stim_start) && (time <= stim_end) && (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
      I_stim = stim_amplitude;
   else
      I_stim = 0.0;

   I_stim = stimulusCurrent;
   
   if (previousSolution[4] > -100.0)
      Xi = 2.837*(exp(0.04*(previousSolution[4]+77.0))-1.0)/((previousSolution[4]+77.0)*exp(0.04*(previousSolution[4]+35.0)));
   else
      Xi = 1.0;

   i_K = g_K*previousSolution[7]*Xi*(previousSolution[4]-E_K);
   alpha_K1 = 1.02/(1.0+exp(0.2385*(previousSolution[4]-E_K1-59.215)));
   beta_K1 = (0.49124*exp(0.08032*(previousSolution[4]+5.476-E_K1))+1.0*exp(0.06175*(previousSolution[4]-(E_K1+594.31))))/(1.0+exp(-0.5143*(previousSolution[4]-E_K1+4.753)));
   K1_infinity = alpha_K1/(alpha_K1+beta_K1);
   i_K1 = g_K1*K1_infinity*(previousSolution[4]-E_K1);
   Kp = 1.0/(1.0+exp((7.488-previousSolution[4])/5.98));
   i_Kp = g_Kp*Kp*(previousSolution[4]-E_Kp);
   computedTimeDerivative[4] = -1.0/C*(I_stim+i_Na+i_si+i_K+i_K1+i_Kp+i_b);
   alpha_d = 0.095*exp(-0.01*(previousSolution[4]-5.0))/(1.0+exp(-0.072*(previousSolution[4]-5.0)));
   beta_d = 0.07*exp(-0.017*(previousSolution[4]+44.0))/(1.0+exp(0.05*(previousSolution[4]+44.0)));
   computedTimeDerivative[5] = alpha_d*(1.0-previousSolution[5])-beta_d*previousSolution[5];
   alpha_f = 0.012*exp(-0.008*(previousSolution[4]+28.0))/(1.0+exp(0.15*(previousSolution[4]+28.0)));
   beta_f = 0.0065*exp(-0.02*(previousSolution[4]+30.0))/(1.0+exp(-0.2*(previousSolution[4]+30.0)));
   computedTimeDerivative[6] = alpha_f*(1.0-previousSolution[6])-beta_f*previousSolution[6];
   alpha_X = 0.0005*exp(0.083*(previousSolution[4]+50.0))/(1.0+exp(0.057*(previousSolution[4]+50.0)));
   beta_X = 0.0013*exp(-0.06*(previousSolution[4]+20.0))/(1.0+exp(-0.04*(previousSolution[4]+20.0)));
   computedTimeDerivative[7] = alpha_X*(1.0-previousSolution[7])-beta_X*previousSolution[7];
}

double ComputeIonicCurrent(double *previousSolution)
{
   if (!cellModelInitialised)
   {
       mexErrMsgIdAndTxt( "MATLAB:luo_rudy:not_initialised",
                "init() must be called before ComputeIonicCurrent");
   }
    
   if (previousSolution==NULL)
   {
       mexErrMsgIdAndTxt( "MATLAB:luo_rudy:null_pointer",
                "previousSolution doesn't point at a valid memory segment");
   }
    
   i_b = g_b*(previousSolution[4]-E_b);
   i_Na = g_Na*pow(previousSolution[2], 3.0)*previousSolution[0]*previousSolution[1]*(previousSolution[4]-E_Na);
   E_si = 7.7-13.0287*log(previousSolution[3]/1.0);
   i_si = 0.09*previousSolution[5]*previousSolution[6]*(previousSolution[4]-E_si);
   if (previousSolution[4] > -100.0)
      Xi = 2.837*(exp(0.04*(previousSolution[4]+77.0))-1.0)/((previousSolution[4]+77.0)*exp(0.04*(previousSolution[4]+35.0)));
   else
      Xi = 1.0;
   i_K = g_K*previousSolution[7]*Xi*(previousSolution[4]-E_K);
   alpha_K1 = 1.02/(1.0+exp(0.2385*(previousSolution[4]-E_K1-59.215)));
   beta_K1 = (0.49124*exp(0.08032*(previousSolution[4]+5.476-E_K1))+1.0*exp(0.06175*(previousSolution[4]-(E_K1+594.31))))/(1.0+exp(-0.5143*(previousSolution[4]-E_K1+4.753)));
   K1_infinity = alpha_K1/(alpha_K1+beta_K1);
   i_K1 = g_K1*K1_infinity*(previousSolution[4]-E_K1);
   Kp = 1.0/(1.0+exp((7.488-previousSolution[4])/5.98));
   i_Kp = g_Kp*Kp*(previousSolution[4]-E_Kp);    
    
   return i_Na+i_si+i_K+i_K1+i_Kp+i_b;
}


//==============================================================================
// End of file
//==============================================================================
