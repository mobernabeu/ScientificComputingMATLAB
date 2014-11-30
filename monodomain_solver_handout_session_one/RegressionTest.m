function RegressionTest( problemConfiguration, conductionVelocities )
%REGRESSIONTEST This function checks the output of a simulation with 50
%cells and conductivities of 1.4 mS/cm against some known conduction
%velocity results

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

recordedConductionVelocity = 0.0521720;

if (problemConfiguration.numberCells==50 && problemConfiguration.leftConductivity==1.4 && problemConfiguration.rightConductivity==1.4)
   if any(conductionVelocities == Inf)
       warning('The simulation was not run for long enough to compute conduction velocities. Regression test could not to be performed.');
       return;
   end
   if any(abs(conductionVelocities - recordedConductionVelocity)>1e-4)
       error('Regression test failed. You may have introduced a bug while optimising the code! Review your last changes.');
   end
end

end

