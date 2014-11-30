function [ ionicCurrents, apCellModelsState ] = CellModelsComputeIonicCurrents( apCellModelsState, problemConfiguration, timeStepNumber, stimuli )
%CELLMODELSCOMPUTEIONICCURRENTS Solves action potential cells models and
%returns transmembrane ionic currents

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

% Work out current time based on time step number
currentTime = problemConfiguration.pdeTimeStep * (timeStepNumber-1);
odeTimeStep = problemConfiguration.pdeTimeStep / problemConfiguration.odeSolvesPerPdeTimeStep;

ionicCurrents = zeros(problemConfiguration.numberCells,1);

odeSolvesPerPdeTimeStep = problemConfiguration.odeSolvesPerPdeTimeStep;

for cellIndex = 1:problemConfiguration.numberCells
    % Compute the ionic current I(apCellModelsState_{n-1}, V_{n-1}). apCellModelsState_{n} means apCellModelsState at time t_{n}
    ionicCurrents(cellIndex,1) = luo_rudy_1991_iionic_mex_adaptor(apCellModelsState(:,cellIndex));
    
    % The evolution of the cell models is given by an ODE of the form 
    % du/dt = f(u,t). We approximate the time derivative with the forward
    % Euler method, i.e. 
    %     d(apCellModelsState)/dt = luo_rudy_1991(apCellModelsState,t)
    %     (apCellModelsState_{n} - apCellModelsState_{n-1})/timeStep ~= luo_rudy_1991(apCellModelsState_{n-1},t_{n})
    %     apCellModelsState_{n} ~= apCellModelsState_{n-1} + timeStep * luo_rudy_1991(apCellModelsState_{n-1},t_{n})
    % The user can choose to perform multiple ODE solves per PDE time step
    % in order to improve accuracy.
    for odeSolveIndex= 1:odeSolvesPerPdeTimeStep
        apCellModelsState(:,cellIndex) = apCellModelsState(:,cellIndex) + (odeTimeStep)*luo_rudy_1991_time_deriv_mex_adaptor(currentTime, ...
                                                                                                                             apCellModelsState(:,cellIndex), ...
                                                                                                                             stimuli(cellIndex));
    end
end

end

