function [ systemRHS ] = AssemblyFiniteDifferencesRightHandSide( problemConfiguration, previousSolution, ionicCurrents, stimuli )
%ASSEMBLYFINITEDIFFERENCESRIGHTHANDSIDE Assembly a semi-implicit finite
%difference monodomain right-hand-side vector according to
%problemConfiguration and previousSolution
%
%  systemRHS = a*b*x - a*c*y + c*z
%
%  Scalar parameters:
%    a = cellSurfaceAreaToVolumeRatio
%    b = cellMembraneCapacitance
%    c = timeStep
%
%  Vector parameters:
%    x = previousSolution
%    y = ionicCurrents
%    z = stimuli

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

systemRHS = zeros(problemConfiguration.numberCells, 1);

for cellIndex=1:problemConfiguration.numberCells
    systemRHS(cellIndex,1) = problemConfiguration.cellSurfaceAreaToVolumeRatio*problemConfiguration.cellMembraneCapacitance*previousSolution(cellIndex) ...
                           - problemConfiguration.cellSurfaceAreaToVolumeRatio*problemConfiguration.pdeTimeStep*ionicCurrents(cellIndex) ...
                           + problemConfiguration.pdeTimeStep*stimuli(cellIndex);
end

end

