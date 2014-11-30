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
    double *previousSolution;
    double ionicCurrent;
     
    /* Check for proper number of arguments */    
    if (nrhs != 1) { 
	    mexErrMsgIdAndTxt( "MATLAB:yprime:invalidNumInputs",
                "One input argument required."); 
    } else if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "MATLAB:yprime:maxlhs",
                "Too many output arguments provided."); 
    } 

    // Do some initialisations required by the cell model.
    init();
    
    // Get a valid C pointer to the first input argument provided
    previousSolution = mxGetPr(prhs[0]);
    
    // Call the C function that computes the ionic current
    ionicCurrent = ComputeIonicCurrent(previousSolution);
    
    // Convert the output into a MATLAB data structure to be returned by the MEX function
    plhs[0] = mxCreateDoubleScalar(ionicCurrent);
    
    return;
}
