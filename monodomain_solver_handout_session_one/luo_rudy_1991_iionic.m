function [ Iionic ] = luo_rudy_1991_iionic( Y )
%LUO_RUDY_1991_IIONIC Summary of this function goes here
%   Detailed explanation goes here

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

E_b = -59.87;   % millivolt (in background_current)
g_b = 0.03921;   % milliS_per_cm2 (in background_current)
i_b = g_b*(Y(5)-E_b);

g_Na = 23.0;   % milliS_per_cm2 (in fast_sodium_current)
Nai = 18.0;   % millimolar (in ionic_concentrations)
Nao = 140.0;   % millimolar (in ionic_concentrations)
F = 96484.6;   % coulomb_per_mole (in membrane)
R = 8314.0;   % joule_per_kilomole_kelvin (in membrane)
T = 310.0;   % kelvin (in membrane)
E_Na = R*T/F*log(Nao/Nai);
i_Na = g_Na*Y(3)^3.0*Y(1)*Y(2)*(Y(5)-E_Na);

E_si = 7.7-13.0287*log(Y(4)/1.0);
i_si = 0.09*Y(6)*Y(7)*(Y(5)-E_si);

Ki = 145.0;   % millimolar (in ionic_concentrations)
Ko = 5.4;   % millimolar (in ionic_concentrations)
g_K1 = 0.6047*sqrt(Ko/5.4);
E_K1 = R*T/F*log(Ko/Ki);
alpha_K1 = 1.02/(1.0+exp(0.2385*(Y(5)-E_K1-59.215)));
beta_K1 = (0.49124*exp(0.08032*(Y(5)+5.476-E_K1))+1.0*exp(0.06175*(Y(5)-(E_K1+594.31))))/(1.0+exp(-0.5143*(Y(5)-E_K1+4.753)));
K1_infinity = alpha_K1/(alpha_K1+beta_K1);
i_K1 = g_K1*K1_infinity*(Y(5)-E_K1);

g_Kp = 0.0183;   % milliS_per_cm2 (in plateau_potassium_current)
Kp = 1.0/(1.0+exp((7.488-Y(5))/5.98));
E_Kp = E_K1;
i_Kp = g_Kp*Kp*(Y(5)-E_Kp);

PR_NaK = 0.01833;   % dimensionless (in time_dependent_potassium_current)
E_K = R*T/F*log((Ko+PR_NaK*Nao)/(Ki+PR_NaK*Nai));
g_K = 0.282*sqrt(Ko/5.4);
if (Y(5) > -100.0)
   Xi = 2.837*(exp(0.04*(Y(5)+77.0))-1.0)/((Y(5)+77.0)*exp(0.04*(Y(5)+35.0)));
else
   Xi = 1.0;
end;
i_K = g_K*Y(8)*Xi*(Y(5)-E_K);

Iionic = i_b + i_Na + i_Kp + i_si + i_K + i_K1;


end

