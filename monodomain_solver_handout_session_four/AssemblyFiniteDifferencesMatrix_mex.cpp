// Introduction to Scientific and High Performance Computing with MATLAB
//
// Copyright (c) 2012-2014, Miguel O. Bernabeu, All rights reserved.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library.

#include "mex.h"
/*
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
*/

#define acces_2d(matrix,rowIndex,columnIndex,numRows) (matrix[(columnIndex)*(numRows) + (rowIndex)])

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{
    double cellSurfaceAreaToVolumeRatio = mxGetScalar( mxGetField(prhs[0], 0, "cellSurfaceAreaToVolumeRatio") );
    double cellMembraneCapacitance = mxGetScalar( mxGetField(prhs[0], 0, "cellMembraneCapacitance") );
    double leftConductivity = mxGetScalar( mxGetField(prhs[0], 0, "leftConductivity") );
    double pdeTimeStep = mxGetScalar( mxGetField(prhs[0], 0, "pdeTimeStep") );
    double spaceStep = mxGetScalar( mxGetField(prhs[0], 0, "spaceStep") );
    double rightConductivity = mxGetScalar( mxGetField(prhs[0], 0, "rightConductivity") );
    unsigned numberCells = mxGetScalar( mxGetField(prhs[0], 0, "numberCells") );
        
    double surfaceRatioTimesCapacitance = cellSurfaceAreaToVolumeRatio * cellMembraneCapacitance;
    double leftCondutivityTimesTimeStepOverSpaceStepSquared = leftConductivity * pdeTimeStep / (spaceStep*spaceStep);
    double rightCondutivityTimesTimeStepOverSpaceStepSquared = rightConductivity * pdeTimeStep / (spaceStep*spaceStep);        
    
    // Create a dense matrix
    mxArray* denseMatrix = mxCreateDoubleMatrix(numberCells, numberCells, mxREAL);
    
    // Get a C pointer to the dense matrix 
    double *systemMatrix = mxGetPr(denseMatrix);

    unsigned rowIndex;
    double condutivityTimesTimeStepOverSpaceStepSquared;
    
    // Fill in matrix row 0
    acces_2d(systemMatrix,0,0,numberCells) = surfaceRatioTimesCapacitance + 2*leftCondutivityTimesTimeStepOverSpaceStepSquared;
    acces_2d(systemMatrix,0,1,numberCells) = -2*leftCondutivityTimesTimeStepOverSpaceStepSquared;

    // Fill in rows 1 to numberCells-2
    for (rowIndex=1;rowIndex<numberCells-1;rowIndex++)
    {
        // Determine which conductivity to use (left half or right half)
        if (rowIndex < (numberCells)/2)
        {
            condutivityTimesTimeStepOverSpaceStepSquared = leftCondutivityTimesTimeStepOverSpaceStepSquared;
        }
        else
        {
            condutivityTimesTimeStepOverSpaceStepSquared = rightCondutivityTimesTimeStepOverSpaceStepSquared;
        }
    
        acces_2d(systemMatrix,rowIndex,rowIndex-1,numberCells) = -condutivityTimesTimeStepOverSpaceStepSquared;
        acces_2d(systemMatrix,rowIndex,rowIndex,numberCells) = surfaceRatioTimesCapacitance + 2*condutivityTimesTimeStepOverSpaceStepSquared;
        acces_2d(systemMatrix,rowIndex,rowIndex+1,numberCells) = -condutivityTimesTimeStepOverSpaceStepSquared;
    }

    // Fill in row numberCells-1
    acces_2d(systemMatrix,numberCells-1,numberCells-2,numberCells) = -2*rightCondutivityTimesTimeStepOverSpaceStepSquared;
    acces_2d(systemMatrix,numberCells-1,numberCells-1,numberCells) = surfaceRatioTimesCapacitance + 2*rightCondutivityTimesTimeStepOverSpaceStepSquared;    
    
    // Do the conversion from dense to sparse (the output goes into plhs, ready to be returned by the MEX function)
    mexCallMATLAB(1,&plhs[0],1,&denseMatrix,"sparse");
    
    // We don't need the dense matrix anymore
    mxDestroyArray(denseMatrix);
    
    return;
}
