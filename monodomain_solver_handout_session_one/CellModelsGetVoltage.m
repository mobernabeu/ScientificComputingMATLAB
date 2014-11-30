function [ voltages ] = CellModelsGetVoltage( problemConfiguration, apCellModelsState )
%CELLMODELSGETVOLTAGE Transmembrane potential is stored as the 5th entry of
%the cell state variable vector. This function gathers the value for all
%the cells in the simulation in a separate vector

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

for cellIndex = 1:problemConfiguration.numberCells
    voltages(cellIndex,1) = apCellModelsState(5,cellIndex);
end

end

