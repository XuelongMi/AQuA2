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
 * ----------------------------------------
 * 
 * This programm demonstrates the usage of the compact watershed implementation
 * from:
 * 
 * "Compact Watershed and Preemptive SLIC:\\On improving trade-offs of superpixel segmentation algorithms"
 * Peer Neubert and Peter Protzel, ICPR 2014 * 
 * 
 */ 

#include "opencv2/opencv.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgproc/imgproc_c.h"
#include "opencv2/core/internal.hpp"
#include <stdio.h>
#include <string>

#include "sys/time.h"
#include "compact_watershed.h"

using namespace std;
using namespace cv;

// start a time measurement
double startTimeMeasure(){
  timeval tim;
  gettimeofday(&tim, 0);
  return tim.tv_sec+(tim.tv_usec/1000000.0);
}

// stop a time measurement
double stopTimeMeasure(double t1){
  timeval tim;
  gettimeofday(&tim, 0);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  return t2-t1;
}

// Overly the image I with a binary mask B and
// store results in R
void getOverlayedImage(Mat& R, Mat& B, Mat& I)
{
  int dimx = I.cols;
  int dimy = I.rows;
  
  // overlayed image  
  I.copyTo(R);
  for(int i=0; i<dimy; i++)
  {
    for(int j=0; j<dimx; j++)
    {
      if( B.at<uchar>(i,j) )
      {
        R.at<cv::Vec3b>(i,j)[0] = 0;
        R.at<cv::Vec3b>(i,j)[1] = 0;
        R.at<cv::Vec3b>(i,j)[2] = 255;
      }
    }
  } 
}

/* Demonstrate the usage of compact watershed on a single image
 * stored on disk.
 * The file is loaded and a superpixel segmentation computed for several 
 * values of the compactness parameter. 
 */
int main(int n, char** vals)
{
  cout << "Compact watershed demo, Peer Neubert, 2014"<<endl;
  cout << "Going to call compact_watershed 300 times"<<endl;
  
  
  // load image
  Mat I = imread("lena.bmp");
    
  for(float compactness=0; compactness<=2; compactness+=1)
  {
    // process
    Mat B_cws, seeds;
    double t=startTimeMeasure();
    for(int i=0; i<100; i++)
      compact_watershed(I, B_cws, 1000, compactness, seeds);
    t = stopTimeMeasure(t) / 100;
    cout<<"CWS: compactness="<<compactness<<",  time in ms: "<<t*1000<<endl;
              
    // overlay image          
    Mat R_cws;
    getOverlayedImage(R_cws, B_cws, I);
    
    // show
    stringstream filename; 
    filename << "cws_compactness_"<<compactness<<".png";
    imshow(filename.str(), R_cws);
  }

  // press any to exit
  int k = waitKey(0);
  
  
}