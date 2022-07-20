/*
 * Compact Watershed
 * Copyright (C) 2014  Peer Neubert, peer.neubert@etit.tu-chemnitz.de
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */ 
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <vector>
#include <iostream>
#include "sys/time.h"

#include "mex_helper.h"
#include "../compact_watershed.h"

using namespace std;
using namespace cv;

/*
 * call: B = mex_compact_watershed(uint8(I), n, compactness, single(seeds))
 *       B = mex_compact_watershed(uint8(I), n, compactness)
 * 
 * I            ... input image 
 * n            ... number of segments
 * compactness  ... compactness parameter, e.g. 1.0 
 * seeds        ... matrix of initial seeds, each row is [i,j], single values (optional)
 *
 * B ... resulting boundary image
 *
 * compile with    
 *  mex mex_compact_watershed.cpp  mex_helper.cpp ../compact_watershed.cpp $(pkg-config --cflags --libs opencv)
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  // ============ parse input ==============
  // create opencv Mat from argument
  cv::Mat I;
  convertMx2Mat(prhs[0], I);
  int n = (int)(mxGetScalar(prhs[1]));          // number of segments
  double compVal = (double)(mxGetScalar(prhs[2])); // compactness parameter
  
  Mat seeds;
  if(nrhs>=4)
    convertMx2Mat(prhs[3], seeds);
    

  // ================== process ================  
  cv::Mat B;
  compact_watershed(I, B, n, compVal, seeds);
  
  // ================ create output ================
  if( nlhs>0)
  {
    convertMat2Mx(B, plhs[0]);  
  }
    
}
  
  
  