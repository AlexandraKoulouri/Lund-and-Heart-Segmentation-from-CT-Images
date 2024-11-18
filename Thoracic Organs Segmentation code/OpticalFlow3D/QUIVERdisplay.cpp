#include <mex.h>
#include <mat.h>
#include <matrix.h>

#define cimg_plugin "cimgmatlab.h"

#include "CImg.h"
#include <iostream>
#include <string>

using namespace cimg_library;
using namespace std;

// The lines below are necessary when using a non-standard compiler as visualcpp6.
#ifdef cimg_use_visualcpp6
#define std
#endif
#ifdef min
#undef min
#undef max
#endif

#ifndef cimg_imagepath
#define cimg_imagepath "img/"
#endif

//global constants

const bool normalize = true; //"Histogram normalization of the images"
const bool morph = true;//"Morphing mode"
const bool imode  = true;//"Complete interpolation (or last frame is missing)"
const bool dispflag = true;//"Visualization"
                

/////////////////////////////////////////////////////////////////////////////////
//Function Prototype
template<typename T> CImg<T> getwarp(const CImg<T>& src, const CImg<>& u);



//mex function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs < 2) mexErrMsgTxt("No enough input arguments.");
  if (nrhs > 2) mexErrMsgTxt("Too many input arguments.");
  if (nrhs == 2){
                CImg<> src(prhs[0],true), u(prhs[1],false);
                          
                //Input images preprocessing
                 CImg<> src_blur(src.size());  
                 src_blur = normalize?src.get_blur(0.5f).equalize(256):src.get_blur(0.5f);
                 
                 CImgDisplay disp(src_blur);
                 
                 if (dispflag) {
                                unsigned int w = src.dimx(), h = src.dimy();
                                const unsigned int dmin = cimg::min(w,h), minsiz = 512;
                                if (dmin<minsiz) { w=w*minsiz/dmin; h=h*minsiz/dmin; }
                                const unsigned int dmax = cimg::max(w,h), maxsiz = 1024;
                                if (dmax>maxsiz) { w=w*maxsiz/dmax; h=h*maxsiz/dmax; }
                                disp.assign(w,h,"Estimated Motion",0);
                             }
                
                 CImg<> c  = src_blur.get_pointwise_norm(1).resize(src.dimx(),src.dimy()).normalize(0,180);
                 c.display(disp);
                // disp.wait();
                 
                 const unsigned char white = 255;
                 c.draw_quiver(u,&white,0.7f,10,-14,0).display(disp);
                 disp.wait();
              
                      
            }
  return;
  
}

// get_warp() : Return the image src warped by the motion field u.
//------------
template<typename T> CImg<T> getwarp(const CImg<T>& src, const CImg<>& u) {
  CImg<T> warp(src);
  cimg_forXY(warp,x,y) warp(x,y) = (T)src.linear_atXY(x - u(x,y,0), y - u(x,y,1));
  return warp;
}
