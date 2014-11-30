function [ apCellModelsState ] = CellModelsInitialise( problemConfiguration )
%CELLMODELSINITIALISE Initialise the cell state variables

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

% This loop initialises the state variables for each of the cell models (taken from luo_rudy_1991.m)
for cellIndex = 1:problemConfiguration.numberCells
    apCellModelsState(:,cellIndex) = [0.982660523699656, 0.989108212766685, 0.00171338077730188, 0.00017948816388306, -84.3801107371, 0.00302126301779861, 0.999967936476325, 0.0417603108167287]';
end

end

