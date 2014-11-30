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

master_lab = 1;
greetings = {'hello', 'hola', 'nihao', 'bonjour'};

spmd
    % If I'm the master lab do one thing
    if (labindex == master_lab) 
        % Receive a greeting from every other worker
        for lab_number=1:numlabs-1
            [data, source, tag] = labReceive('any');
            fprintf('Lab number %d says %s.\n', source, data);
        end
        
        % Return the greeting to every other worker
        labBroadcast(master_lab, greetings{labindex});       
    
    % If I'm not the master lab do something different
    else
        % Send a greeting to the master lab
        labSend(greetings{labindex}, master_lab);
        
        % Receive a greeting from the master lab
        data = labBroadcast(master_lab);
        fprintf('Lab number %d replies %s.\n', master_lab, data);
    
    end
end
