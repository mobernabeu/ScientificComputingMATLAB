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
#include "luo_rudy_1991_c_implementation.h"

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{
    double time;
    double *previousSolution;
    double stimulusCurrent;
    mwSize m, n;
    double *computedTimeDerivative;

    // Check for proper number of arguments
    if (nrhs != 3) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidNumInputs",
                "Three input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "MATLAB:yprime:maxlhs",
                "Too many output arguments requested."); 
    } 

    // Do some initialisations required by the cell model.
    init();
   
    /* 
     * Use the function mxGetScalar to get the value of the first input argument provided. 
     * Remember that C arrays are indexed from 0.
     */
    time = mxGetScalar(prhs[0]);

    /* 
     * Use the function mxGetPr to get a pointer to the second input argument provided.
     * Remember that C arrays are indexed from 0.
     */
    previousSolution = mxGetPr(prhs[1]);

    /* 
     * Use the function mxGetScalar to get the value of the third input argument provided. 
     * Remember that C arrays are indexed from 0.
     */
    stimulusCurrent = mxGetScalar(prhs[2]);    
    
    // Use mxGetM and mxGetN to obtain the dimensions of the second input argument.
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    
    /*
     * Create a matrix of size (m,n) and type mxREAL with mxCreateDoubleMatrix. This will
     * be the output of the mexFunction.
     */
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    
    // Obtain a valid C pointer to plhs[0]
    computedTimeDerivative = mxGetPr(plhs[0]);

    // Call the function that computes the time derivative
    ComputeTimeDerivative(time,previousSolution,computedTimeDerivative,stimulusCurrent);    
    
    return;
}
