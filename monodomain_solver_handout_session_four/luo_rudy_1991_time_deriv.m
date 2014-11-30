function [ dY ] = luo_rudy_1991_time_deriv(time, Y, iStim)

% Introduction to Scientific and High Performance Computing with MATLAB
%
% Copyright (c) 2012-2014, Miguel O. Bernabeu, All rights reserved.
% 
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 3.0 of the License, or (at your option) any later version.
% 
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public
% License along with this library.

%===============================================================================
% CellML file:   C:\Users\Miguel Bernabeu\Desktop\luo_rudy_1991.cellml
% CellML model:  luo_rudy_1991
% Date and time: 4/13/2012 at 11:30:06 at AM
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2012 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================

%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [0.982660523699656, 0.989108212766685, 0.00171338077730188, 0.00017948816388306, -84.3801107371, 0.00302126301779861, 0.999967936476325, 0.0417603108167287];

% YNames = {'h', 'j', 'm', 'Cai', 'V', 'd', 'f', 'X'};
% YUnits = {'dimensionless', 'dimensionless', 'dimensionless', 'millimolar', 'millivolt', 'dimensionless', 'dimensionless', 'dimensionless'};
% YComponents = {'fast_sodium_current_h_gate', 'fast_sodium_current_j_gate', 'fast_sodium_current_m_gate', 'intracellular_calcium_concentration', 'membrane', 'slow_inward_current_d_gate', 'slow_inward_current_f_gate', 'time_dependent_potassium_current_X_gate'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: h (dimensionless) (in fast_sodium_current_h_gate)
% 2: j (dimensionless) (in fast_sodium_current_j_gate)
% 3: m (dimensionless) (in fast_sodium_current_m_gate)
% 4: Cai (millimolar) (in intracellular_calcium_concentration)
% 5: V (millivolt) (in membrane)
% 6: d (dimensionless) (in slow_inward_current_d_gate)
% 7: f (dimensionless) (in slow_inward_current_f_gate)
% 8: X (dimensionless) (in time_dependent_potassium_current_X_gate)

%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

E_b = -59.87;   % millivolt (in background_current)
g_b = 0.03921;   % milliS_per_cm2 (in background_current)
g_Na = 23.0;   % milliS_per_cm2 (in fast_sodium_current)
Ki = 145.0;   % millimolar (in ionic_concentrations)
Ko = 5.4;   % millimolar (in ionic_concentrations)
Nai = 18.0;   % millimolar (in ionic_concentrations)
Nao = 140.0;   % millimolar (in ionic_concentrations)
C = 1.0;   % microF_per_cm2 (in membrane)
F = 96484.6;   % coulomb_per_mole (in membrane)
R = 8314.0;   % joule_per_kilomole_kelvin (in membrane)
T = 310.0;   % kelvin (in membrane)
stim_amplitude = -25.5;   % microA_per_cm2 (in membrane)
stim_duration = 2.0;   % millisecond (in membrane)
stim_end = 9000.0;   % millisecond (in membrane)
stim_period = 1000.0;   % millisecond (in membrane)
stim_start = 100.0;   % millisecond (in membrane)
g_Kp = 0.0183;   % milliS_per_cm2 (in plateau_potassium_current)
PR_NaK = 0.01833;   % dimensionless (in time_dependent_potassium_current)

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% i_b (microA_per_cm2) (in background_current)
% alpha_h (per_millisecond) (in fast_sodium_current_h_gate)
% beta_h (per_millisecond) (in fast_sodium_current_h_gate)
% alpha_j (per_millisecond) (in fast_sodium_current_j_gate)
% beta_j (per_millisecond) (in fast_sodium_current_j_gate)
% alpha_m (per_millisecond) (in fast_sodium_current_m_gate)
% beta_m (per_millisecond) (in fast_sodium_current_m_gate)
% E_Na (millivolt) (in fast_sodium_current)
% i_Na (microA_per_cm2) (in fast_sodium_current)
% I_stim (microA_per_cm2) (in membrane)
% E_Kp (millivolt) (in plateau_potassium_current)
% Kp (dimensionless) (in plateau_potassium_current)
% i_Kp (microA_per_cm2) (in plateau_potassium_current)
% alpha_d (per_millisecond) (in slow_inward_current_d_gate)
% beta_d (per_millisecond) (in slow_inward_current_d_gate)
% alpha_f (per_millisecond) (in slow_inward_current_f_gate)
% beta_f (per_millisecond) (in slow_inward_current_f_gate)
% E_si (millivolt) (in slow_inward_current)
% i_si (microA_per_cm2) (in slow_inward_current)
% alpha_X (per_millisecond) (in time_dependent_potassium_current_X_gate)
% beta_X (per_millisecond) (in time_dependent_potassium_current_X_gate)
% Xi (dimensionless) (in time_dependent_potassium_current_Xi_gate)
% E_K (millivolt) (in time_dependent_potassium_current)
% g_K (milliS_per_cm2) (in time_dependent_potassium_current)
% i_K (microA_per_cm2) (in time_dependent_potassium_current)
% K1_infinity (dimensionless) (in time_independent_potassium_current_K1_gate)
% alpha_K1 (per_millisecond) (in time_independent_potassium_current_K1_gate)
% beta_K1 (per_millisecond) (in time_independent_potassium_current_K1_gate)
% E_K1 (millivolt) (in time_independent_potassium_current)
% g_K1 (milliS_per_cm2) (in time_independent_potassium_current)
% i_K1 (microA_per_cm2) (in time_independent_potassium_current)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (millisecond)

i_b = g_b*(Y(5)-E_b);
E_Na = R*T/F*log(Nao/Nai);
i_Na = g_Na*Y(3)^3.0*Y(1)*Y(2)*(Y(5)-E_Na);

if (Y(5) < -40.0)
   alpha_h = 0.135*exp((80.0+Y(5))/-6.8);
else
   alpha_h = 0.0;
end;

if (Y(5) < -40.0)
   beta_h = 3.56*exp(0.079*Y(5))+310000.0*exp(0.35*Y(5));
else
   beta_h = 1.0/(0.13*(1.0+exp((Y(5)+10.66)/-11.1)));
end;

dY(1, 1) = alpha_h*(1.0-Y(1))-beta_h*Y(1);

if (Y(5) < -40.0)
   alpha_j = (-127140.0*exp(0.2444*Y(5))-0.00003474*exp(-0.04391*Y(5)))*(Y(5)+37.78)/(1.0+exp(0.311*(Y(5)+79.23)));
else
   alpha_j = 0.0;
end;

if (Y(5) < -40.0)
   beta_j = 0.1212*exp(-0.01052*Y(5))/(1.0+exp(-0.1378*(Y(5)+40.14)));
else
   beta_j = 0.3*exp(-0.0000002535*Y(5))/(1.0+exp(-0.1*(Y(5)+32.0)));
end;

dY(2, 1) = alpha_j*(1.0-Y(2))-beta_j*Y(2);
alpha_m = 0.32*(Y(5)+47.13)/(1.0-exp(-0.1*(Y(5)+47.13)));
beta_m = 0.08*exp(-Y(5)/11.0);
dY(3, 1) = alpha_m*(1.0-Y(3))-beta_m*Y(3);
E_si = 7.7-13.0287*log(Y(4)/1.0);
i_si = 0.09*Y(6)*Y(7)*(Y(5)-E_si);
dY(4, 1) = -0.0001/1.0*i_si+0.07*(0.0001-Y(4));

if ((time >= stim_start) && (time <= stim_end) && (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
   I_stim = stim_amplitude;
else
   I_stim = 0.0;
end;

% Overwrite the built-in stimulus with the one provided by the user
I_stim = iStim;

g_K = 0.282*sqrt(Ko/5.4);

if (Y(5) > -100.0)
   Xi = 2.837*(exp(0.04*(Y(5)+77.0))-1.0)/((Y(5)+77.0)*exp(0.04*(Y(5)+35.0)));
else
   Xi = 1.0;
end;

E_K = R*T/F*log((Ko+PR_NaK*Nao)/(Ki+PR_NaK*Nai));
i_K = g_K*Y(8)*Xi*(Y(5)-E_K);
g_K1 = 0.6047*sqrt(Ko/5.4);
E_K1 = R*T/F*log(Ko/Ki);
alpha_K1 = 1.02/(1.0+exp(0.2385*(Y(5)-E_K1-59.215)));
beta_K1 = (0.49124*exp(0.08032*(Y(5)+5.476-E_K1))+1.0*exp(0.06175*(Y(5)-(E_K1+594.31))))/(1.0+exp(-0.5143*(Y(5)-E_K1+4.753)));
K1_infinity = alpha_K1/(alpha_K1+beta_K1);
i_K1 = g_K1*K1_infinity*(Y(5)-E_K1);
Kp = 1.0/(1.0+exp((7.488-Y(5))/5.98));
E_Kp = E_K1;
i_Kp = g_Kp*Kp*(Y(5)-E_Kp);
dY(5, 1) = -1.0/C*(I_stim+i_Na+i_si+i_K+i_K1+i_Kp+i_b);
alpha_d = 0.095*exp(-0.01*(Y(5)-5.0))/(1.0+exp(-0.072*(Y(5)-5.0)));
beta_d = 0.07*exp(-0.017*(Y(5)+44.0))/(1.0+exp(0.05*(Y(5)+44.0)));
dY(6, 1) = alpha_d*(1.0-Y(6))-beta_d*Y(6);
alpha_f = 0.012*exp(-0.008*(Y(5)+28.0))/(1.0+exp(0.15*(Y(5)+28.0)));
beta_f = 0.0065*exp(-0.02*(Y(5)+30.0))/(1.0+exp(-0.2*(Y(5)+30.0)));
dY(7, 1) = alpha_f*(1.0-Y(7))-beta_f*Y(7);
alpha_X = 0.0005*exp(0.083*(Y(5)+50.0))/(1.0+exp(0.057*(Y(5)+50.0)));
beta_X = 0.0013*exp(-0.06*(Y(5)+20.0))/(1.0+exp(-0.04*(Y(5)+20.0)));
dY(8, 1) = alpha_X*(1.0-Y(8))-beta_X*Y(8);

%===============================================================================
% End of file
%===============================================================================
