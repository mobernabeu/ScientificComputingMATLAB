function my_simulation( input_args )
%MY_SIMULATION Example simulation for course day 1
%   It illustrates the use of tic/toc to time MATLAB code 
%   and more advanced profiling techniques.

% my_simulation has a setup stage and a run stage to be executed in order
tic;
setup_my_simulation();
toc;

tic;
profile on
run_my_simulation();
profile off
toc;


end

