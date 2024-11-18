/*-----------------------------------------------------------------------
# This is a modified version of the following file:
#  
#   File        : optical_flow.cpp 
#   Description : Compute the optical flow between two images, with a multiscale and variational algorithm  
#   Author      : Tschumperlé David 
#   Institution : ODYSSEE, INRIA Sophia Antipolis. 
#   Contact     : David.Tschumperle@sophia.inria.fr 
#   Date        : 03/12/2003 
#    
#   This program is free software; you can redistribute it and/or modify 
#   it under the terms of the GNU General Public License as published by 
#   the Free Software Foundation; either version 2 of the License, or 
#   (at your option) any later version. 
//////////////////////////////////////////////////////////////////////////////////////////
#The Input parameters are two images (2D or 3D) and as an output we have the displacement vector u.
--------------------------------------------------------------------------*/

#include <mex.h>
#include <mat.h>
#include <matrix.h>

#define cimg_plugin "cimgmatlab.h"

#include "CImg.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace cimg_library;
using namespace std;


//global constants
const float smooth = 0.1f;//"Flow Smoothness"
const float precision = 0.09f;//"Convergence precision"
const unsigned int nb=110;//"Number of warped frames"
const unsigned int nbscale = 0 ;//"Number of scales (0=auto)");

const bool normalize = true; //"Histogram normalization of the images"
const bool morph = true;//"Morphing mode"
const bool imode  = true;//"Complete interpolation (or last frame is missing)"
const bool dispflag = true;//"Visualization"
//Functions
CImg<> optmonoflow(const CImg<>& I1, const CImg<>& I2, const CImg<>& u0,
                   const float smooth, const float precision, CImgDisplay& disp);
CImg<> optflow(const CImg<>& xsrc, const CImg<>& xdest,
               const float smooth, const float precision, const unsigned int pnb_scale, CImgDisplay& disp);
                
//mex function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs < 2) mexErrMsgTxt("No enough input arguments.");
  if (nrhs > 2) mexErrMsgTxt("Too many input arguments.");
  if (nrhs == 2){
                CImg<> src_blur(prhs[0],false), dest_blur(prhs[1],false);//2 Input Volumes
                //Input images preprocessing
                // src_blur = normalize?src_blur.get_blur(0.5f).equalize(256):src_blur.get_blur(0.5f),
               //  dest_blur  = normalize?dest_blur.get_blur(0.5f).equalize(256):dest_blur.get_blur(0.5f);
                 CImgDisplay disp;
                 const CImg<> u = optflow(src_blur,dest_blur,smooth,precision,nbscale,disp);
                
                 plhs[0] = u.toMatlab();
                 ofstream outputFilex, outputFiley,outputFilez;
                 outputFilex.open("displacement_x.txt");
                 outputFiley.open("displacement_y.txt");
                 outputFilez.open("displacement_z.txt");
                 cimg_forY(u,y){
                      cimg_forXZ(u,x,z){
                    //float mag = (float) sqrt(pow(u(x,y,z,0),2)+pow(u(x,y,z,1),2)+pow(u(x,y,z,2),2));
                          outputFilex<<u(x,y,z,0)<<" ";
                          outputFiley<<u(x,y,z,1)<<" ";
                          outputFilez<<u(x,y,z,2)<<" ";
                     }
                  outputFilex<<endl;outputFiley<<endl;outputFilez<<endl;
                 }
                outputFilex.close();outputFiley.close();outputFilez.close();
            
 }
  return;
     
}
/////////////////////////////////////////////////////////////////////////////////

CImg<> optflow(const CImg<>& xsrc, const CImg<>& xdest,
               const float smooth, const float precision, const unsigned int pnb_scale, CImgDisplay& disp) {
  const CImg<>
    src  = xsrc.get_pointwise_norm(1).resize(xdest.dimx(),xdest.dimy(),xdest.dimz()).normalize(0,1),
    dest = xdest.get_pointwise_norm(1).resize(xdest.dimx(),xdest.dimy(),xdest.dimz()).normalize(0,1);
    CImg<> u =  CImg<>(src.dimx(),src.dimy(),src.dimz(),3,1).fill(0);
       
  //multi-scalling 
    const unsigned int nb_scale = pnb_scale>0?pnb_scale:(unsigned int)(2*log((double)(max(max(src.dimx(),src.dimy()),src.dimz()))));
  //const unsigned int nb_scale = pnb_scale>0?pnb_scale:(unsigned int)(2*std::log((double)(cimg::max(src.dimx(),src.dimy()))));
  for (int scale=nb_scale-1; scale>=0; scale--) {
       const CImg<> I1 = src.get_resize((int)(ceil(src.dimx()/std::pow(1.5,scale))), (int)(ceil(src.dimy()/std::pow(1.5,scale))) ,(int)(ceil(src.dimz()/std::pow(1.5,scale))));
       const CImg<> I2 = dest.get_resize((int)(ceil(src.dimx()/std::pow(1.5,scale))), (int)(ceil(src.dimy()/std::pow(1.5,scale))) ,(int)(ceil(src.dimz()/std::pow(1.5,scale))));
       u*=1.5;
       u = optmonoflow(I1,I2,u,smooth,(float)(precision/std::pow(2.25,1+scale)),disp);//apply optical flow algorithm for different scales
    
  }
    cimg_forXYZV(u,x,y,z,v){
    if(src(x,y,z)==0) u(x,y,z,v)=0;}
 return u;
}

// optmonoflow() : estimate optical flow for one scale ( semi-implicite PDE scheme ) between I2->I1
//---------------
CImg<> optmonoflow(const CImg<>& I1, const CImg<>& I2, const CImg<>& u0,
                   const float smooth, const float precision, CImgDisplay& disp) {

 CImg<float> u = u0.get_resize(I1.dimx(),I1.dimy(),I1.dimz(),3);
 CImg<float> dI(I2.dimx(),I2.dimy(),I2.dimz(),3);
  dI.fill(0);
 //CImg_3x3x3(I,float);
  float dt=2,E=1e20f;

  // compute first derivatives of I2
    cimg_for3XYZ(dI,x,y,z){
              dI(x,y,z,0) = 0.5*(I2(_n1x,y,z)-I2(_p1x,y,z));//derivatives of I2
              dI(x,y,z,1) = 0.5*(I2(x,_n1y,z)-I2(x,_p1y,z));
              dI(x,y,z,2) = 0.5*(I2(x,y,_n1z)-I2(x,y,_p1z));
                          
     }
   
// Main PDE iteration
for (unsigned int iter=0; iter<50; iter++) {
    std::fprintf(stderr,"\r- Iteration %d - E = %g",iter,E); std::fflush(stderr);
    const float Eold = E;
    E = 0;
    cimg_for3XYZ(u,x,y,z) { //3x3 neighborhood  
      const float
        X = x + u(x,y,z,0),
        Y = y + u(x,y,z,1),
        Z = z + u(x,y,z,2),
        deltaI = (float)(I2.linear_atXYZ(X,Y,Z) - I1(x,y,z));// partial derivative over time-dI/dt
      float tmpf = 0;
      cimg_forV(u,k) {
        const float
          ux  = 0.5f*(u(_n1x,y,z,k)-u(_p1x,y,z,k)),//derivatives of velocity(displacemenet)
          uy  = 0.5f*(u(x,_n1y,z,k)-u(x,_p1y,z,k)),
          uz  = 0.5f*(u(x,y,_n1z,k)-u(x,y,_p1z,k)); 
          u(x,y,z,k) = (float)( u(x,y,z,k) +
                            dt*(
                                -deltaI*dI.linear_atXYZ(X,Y,Z,k) +
                                smooth* ( u(_n1x,y,z,k) + u(_p1x,y,z,k) + u(x,_n1y,z,k) + u(x,_p1y,z,k)+  u(x,y,_n1z,k) + u(x,y,_p1z,k)-6*u(x,y,z,k))
                                )/(1+4*smooth*dt)
                              );
      tmpf += ux*ux + uy*uy + uz*uz;
     
      }
      E += deltaI*deltaI + smooth * tmpf;
    }
    if (cimg::abs(Eold-E)<precision) break;
    if (Eold<E) dt*=0.5;
 
}
  return  u;
}