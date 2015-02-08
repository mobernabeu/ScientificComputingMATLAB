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

%
% This sript runs two monodomain simulations with different input
% parameters on MATLAB's default parallel job manager
%

% Creates a list of files that need to be copied over when running in
% remote resources
required_files = { 'AssemblyFiniteDifferencesMatrix.m' 
                   'AssemblyFiniteDifferencesMatrix_mex.mexa64' 
                   'AssemblyFiniteDifferencesRightHandSide.m' 
                   'CellModelsComputeIonicCurrents.m' 
                   'CellModelsGetVoltage.m' 
                   'CellModelsInitialise.m' 
                   'CellModelsSetVoltage.m' 
                   'GetStimuliForTimeStep.m' 
                   'RegressionTest.m' 
                   'RunAndVisualiseMonodomainSimulation.m' 
                   'SolveLinearSystem.m' 
                   'luo_rudy_1991_iionic.m' 
                   'luo_rudy_1991_iionic_mex_adaptor.mexa64' 
                   'luo_rudy_1991_time_deriv.m'
                   'luo_rudy_1991_time_deriv_mex_adaptor.mexa64'};

% Obtains a reference to the scheduler defined on MATLAB's default parallel
% configuration.
sched = parcluster();

% Obtains a reference to the Legion scheduler
%sched = parcluster('legion_R2013a');

% Creates a job that will be submitted to the selected scheduler.
job = createJob(sched);

% Set the files that need to be copied to the cluster backend.
%job.AttachedFiles = required_files;

% Creates a MATLAB task to be executed as part of the job. It will consist
% of a run of function RunAndVisualiseMonodomainSimulation. The rest of
% arguments indicate that it returns two parameters and takes a cell vector
% with the input parameters 
task = createTask(job, @RunAndVisualiseMonodomainSimulation,3,{50,1.5,1.4,1.4,false});

% Submit the job
submit(job);

% Create and submit a new job with a different set of input parameters
job2 = createJob(sched);
%job2.AttachedFiles = required_files;
createTask(job2, @RunAndVisualiseMonodomainSimulation,3,{120,1.5,1.4,1.4,false});
submit(job2);

% Wait for both jobs to finish
wait(job,'finished');
wait(job2,'finished');

% Copy the results over from the remote resource
job_results  = fetchOutputs(job);
job2_results = fetchOutputs(job2);

% Display the results
celldisp(job_results);
celldisp(job2_results);

% Do some postprocessing
% EDIT HERE

% Once the results have been retrieved we can remove the job from the queue
destroy(job);
destroy(job2);
