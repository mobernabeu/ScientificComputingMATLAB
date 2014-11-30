function [ stimuli ] = GetStimuliForTimeStep( problemConfiguration, timeStepNumber )
%GETSTIMULIFORTIMESTEP Returns the volumetric stimulus to be applied at
%time step timeStepNumber

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
currentTime = problemConfiguration.pdeTimeStep*(timeStepNumber-1);

% Initialise stimuli vector to zero
stimuli = zeros(problemConfiguration.numberCells, 1);

% If its time to do it, stimulate two cells in the middle of the domain
if (currentTime>=problemConfiguration.stimulusStart && currentTime<problemConfiguration.stimulusStart+problemConfiguration.stimulusDuration)
    stimuli(floor(problemConfiguration.numberCells/2) : floor(problemConfiguration.numberCells/2) + 1) = problemConfiguration.stimulusMagnitude;
end

end

