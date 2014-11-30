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

spmd
    % We create a codistributed array representing an identity matrix
    eye_matrix = codistributed.eye(10);    
    
    % We compute the eigenvalues of the identity matrix (all ones).
    % The result is also a codistributed array
    eigenvalues = eig(eye_matrix);
    
    % This variable is local to each worker.
    num_eigenvalues_per_worker = size(getLocalPart(eigenvalues),1);
end

% After the spmd region, eigenvalues becomes a distributed array and the
% client has complete access to it.
eigenvalues

for lab_id=1:matlabpool('size')
    % num_eigenvalues_per_worker has become a composite object after
    % leaving the spmd region.
    num_eig = num_eigenvalues_per_worker{lab_id};
    fprintf('Lab number %d stored %d eigenvalues.\n', lab_id, num_eig);
end
