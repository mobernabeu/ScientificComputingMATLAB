function [ finalVoltage, conductionVelocities ] = RunAndVisualiseMonodomainSimulation( numberCells, simulationDuration, leftHalfFibreConductivity, rightHalfFibreConductivity, visualiseSimulation)
%RUNANDVISUALISEMONODOMAINSIMULATION runs, and optionally visualises, a
%monodomain simulation in a one-dimensional domain (cardiac tissue fibre).
%
%  [finalVoltage, conductionVelocities] =
%                 RUNANDVISUALISEMONODOMAINSIMULATION(numberCells,
%                                                     simulationDuration,
%                                                     leftFibreHalfconductivity,
%                                                     rightFibreHalfconductivity,
%                                                     visualiseSimulation)
%
%  Runs a monodomain simulation in a cardiac fibre made of numberCells
%  cells for simulationDuration milliseconds. A square stimulus is applied
%  after 0.5 ms of the start of the simulation for 1 ms at the middle of
%  the domain. leftHalfFibreConductivity and rightHalfFibreConductivity are
%  the conductivities used in the left and right half of the fibre,
%  respectively. If visualiseSimulation is true, the evolution of the
%  transmembrane potential throughout the simulation is plotted, as well as
%  the action potential shape at both ends of the fibre.
%
%  The simulation returns the transmembrane potential accross the domain
%  at the end of the simulation (finalVoltage) and a vector
%  with the conduction velocity computed at each end of the fibre
%  (conductionVelocities).

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
% Section 1: setup simulation
%
problemConfiguration.spaceStep = 0.025; % distance between cells (cm)
problemConfiguration.pdeTimeStep = 0.01; % monodomain temporal discretisation time step (ms)
problemConfiguration.odeSolvesPerPdeTimeStep = 10; % number of ODE solves per pde time step (i.e. odeTimeStep = pdeTimeStep/odeSolvesPerPdeTimeStep)
problemConfiguration.stimulusMagnitude = 1e5; % magnitude of the external stimulus applies (A/cm^3)
problemConfiguration.stimulusStart = 0.5; % stimulus application delay (ms)
problemConfiguration.stimulusDuration = 1; % stimulus duration (ms)
problemConfiguration.leftConductivity = leftHalfFibreConductivity; % conductivity in the left half of the fibre (mS/cm)
problemConfiguration.rightConductivity = rightHalfFibreConductivity; % conductivity in the right half of the fibre (mS/cm)
problemConfiguration.cellSurfaceAreaToVolumeRatio = 1400; % cell surface-area to volume ratio (cm^-1)
problemConfiguration.cellMembraneCapacitance = 1; % cell membrane capacitance (uF/cm^2)
problemConfiguration.numberCells = numberCells; % number of cells in the domain

apCellModelsState = CellModelsInitialise(problemConfiguration); % initialise matrix storing column-wise the state of each action potential cell model

previousSolution = CellModelsGetVoltage(problemConfiguration, apCellModelsState);  % Get resting potential from cell models

numberTimeSteps = simulationDuration / problemConfiguration.pdeTimeStep; % Number of time steps to be executed
%
% End Section 1
%

%
% Section 6.a: visualise results
%
if visualiseSimulation
    gridPoints = (1:numberCells)';
    figure(1);
    plot(gridPoints,previousSolution,'b');
    xlabel('Cell number'); ylabel('Transmembrane potential (mV)'); title('Transmembrane potential along the fibre, timestep 0');
    axis([1 numberCells -90.0 70.0]);
    drawnow; % plot initial solution
end

%
% Section 2: assembly monodomain system matrix
%
systemMatrix = AssemblyFiniteDifferencesMatrix(problemConfiguration);
% End Section 2

firstCellActivationTime = 0; lastCellActivationTime = 0;
firstCellActive = false;lastCellActive = false;
for timeStepNumber = 1:numberTimeSteps
    %
    % Section 3: solve cell models to compute ionic currents
    %
    volumetricStimuli = GetStimuliForTimeStep(problemConfiguration, timeStepNumber);
    [ionicCurrents, apCellModelsState] = CellModelsComputeIonicCurrents(apCellModelsState, ...
                                                                        problemConfiguration, ...
                                                                        timeStepNumber, ...
                                                                        volumetricStimuli/problemConfiguration.cellSurfaceAreaToVolumeRatio);
    
    %
    % Section 4: assembly monodomain system right hand side
    %
    systemRightHandSide = AssemblyFiniteDifferencesRightHandSide(problemConfiguration, previousSolution, ionicCurrents, volumetricStimuli);
    
    %
    % Section 5: solve monodomain system
    %
    currentSolution = SolveLinearSystem(systemMatrix, systemRightHandSide, previousSolution);
    apCellModelsState = CellModelsSetVoltage(apCellModelsState, problemConfiguration, currentSolution); % Feed into the cell models the voltage computed at tissue level
    
    %
    % Section 6.b: visualise results
    %
    if visualiseSimulation && mod(timeStepNumber,25)==0
        plot(gridPoints,currentSolution,'b');
        xlabel('Cell number'); ylabel('Transmembrane potential (mV)'); title(['Transmembrane potential along the fibre, timestep ' num2str(timeStepNumber)]);
        axis([1 numberCells -90.0 70.0]);
        drawnow;
    end
    
    %
    % Section 7: postprocess results
    %
    firstCellVoltage(timeStepNumber) = currentSolution(1);
    lastCellVoltage(timeStepNumber) = currentSolution(end);
    if currentSolution(1) > 0 && not(firstCellActive)
        firstCellActivationTime = (timeStepNumber-1) * problemConfiguration.pdeTimeStep;
        firstCellActive = true;
    end
    if currentSolution(end) > 0 && not(lastCellActive)
        lastCellActivationTime = (timeStepNumber-1) * problemConfiguration.pdeTimeStep;
        lastCellActive = true;
    end
    %
    % End Section 7
    %
    
    previousSolution = currentSolution;
end

%
% Section 6.c: visualise results
%
if visualiseSimulation && numberTimeSteps > 0
    figure(2);
    plot(1:numberTimeSteps, firstCellVoltage);
    xlabel('Time step number'); ylabel('Transmembrane potential (mV)'); title('Transmembrane potential over time for cell number 1');
    axis([1 numberTimeSteps -90.0 50.0]);
    drawnow;
    
    figure(3);
    plot(1:numberTimeSteps, lastCellVoltage);
    xlabel('Time step number'); ylabel('Transmembrane potential (mV)'); title(['Transmembrane potential over time for cell number ' num2str(problemConfiguration.numberCells)]);
    axis([1 numberTimeSteps -90.0 50.0]);
    drawnow;
end


halfTheFibreLength = (numberCells - 1) * problemConfiguration.spaceStep / 2;
conductionVelocities = halfTheFibreLength * [1/firstCellActivationTime, 1/lastCellActivationTime];
finalVoltage = previousSolution;

%
% Run a regression test to make sure that bugs have not been introduced
% during optimisation
%
RegressionTest(problemConfiguration, conductionVelocities);

end

