function [ systemMatrix ] = AssemblyFiniteDifferencesMatrix( problemConfiguration )
%ASSEMBLYFINITEDIFFERENCESMATRIX Assembly a semi-implicit finite difference
%monodomain matrix (systemMatrix) according to the problem description in
%problemConfiguration.
%
%   systemMatrix is a tridiagonal matrix of the form
%
%     [a 2*b 0 ...       0]
%     [b  a  b  0 ...    0]
%     [0  b  a  b  0 ... 0]
%     [0 ... b  a  b ... 0]
%     [0 ... 0  b  a  b  0]
%     [0    ... 0  b  a  b]
%     [0       ... 0 2*b a]
%
%   with a = surfaceRatio*capacitance + 2*conductivity*timeStep/spaceStep^2,
%        b = -conductivity*timeStep/spaceStep^2
%
%        conductivity = { leftConductivity,  if columnIndex <= (numberCells+1)/2,
%                       { rightConductivity, otherwise.

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

surfaceRatioTimesCapacitance = problemConfiguration.cellSurfaceAreaToVolumeRatio * problemConfiguration.cellMembraneCapacitance;
leftCondutivityTimesTimeStepOverSpaceStepSquared = problemConfiguration.leftConductivity * problemConfiguration.pdeTimeStep / problemConfiguration.spaceStep^2;
rightCondutivityTimesTimeStepOverSpaceStepSquared = problemConfiguration.rightConductivity * problemConfiguration.pdeTimeStep / problemConfiguration.spaceStep^2;

assert(problemConfiguration.numberCells > 1, 'Minimum problem size is currently 2');

% First matrix row
systemMatrix(1,1) = surfaceRatioTimesCapacitance + 2*leftCondutivityTimesTimeStepOverSpaceStepSquared;
systemMatrix(1,2) = -2*leftCondutivityTimesTimeStepOverSpaceStepSquared;

% Rows 2 to numberCells-1
for rowIndex = 2:problemConfiguration.numberCells-1
    % Determine which conductivity to use (left half or right half)
    if rowIndex <= (problemConfiguration.numberCells+1)/2
        condutivityTimesTimeStepOverSpaceStepSquared = leftCondutivityTimesTimeStepOverSpaceStepSquared;
    else
        condutivityTimesTimeStepOverSpaceStepSquared = rightCondutivityTimesTimeStepOverSpaceStepSquared;
    end
    
    systemMatrix(rowIndex,rowIndex-1) = -condutivityTimesTimeStepOverSpaceStepSquared;
    systemMatrix(rowIndex,rowIndex) = surfaceRatioTimesCapacitance + 2*condutivityTimesTimeStepOverSpaceStepSquared;
    systemMatrix(rowIndex,rowIndex+1) = -condutivityTimesTimeStepOverSpaceStepSquared;
end

% Last matrix row
systemMatrix(problemConfiguration.numberCells,problemConfiguration.numberCells-1) = -2*rightCondutivityTimesTimeStepOverSpaceStepSquared;
systemMatrix(problemConfiguration.numberCells,problemConfiguration.numberCells) = surfaceRatioTimesCapacitance + 2*rightCondutivityTimesTimeStepOverSpaceStepSquared;

end
