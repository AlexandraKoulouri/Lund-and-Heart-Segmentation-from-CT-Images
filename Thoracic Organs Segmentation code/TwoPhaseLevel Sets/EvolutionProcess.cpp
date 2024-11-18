/*/////////////////////////////////////////////////////////////////////////////////////
                            2 phase level set method
 0-input image
 1-initial phi1
 2-dt=>time-step of iteration
 3-kappa =>coefficient of the weighted length term L(phi1)
 4-lambda1 =>coefficient of insice C term
 5-lambda2 =>coefficient of outside C term
 6-lambda=>coefficient of the weighted length term L(phi2)
 7-mu => coefficient of the internal (penalizing) energy term P(phi1),P(phi2)
 8-v=>coefficient of the weighted area term A(phi2)
 9-iterations
 10-g=>edge indicator
/*/////////////////////////////////////////////////////////////////////////////////////
#include <mex.h>
#include <mat.h>
#include <matrix.h>

#define cimg_plugin "cimgmatlab.h"

#include "CImg.h"
#include <iostream>
#include <string>
#include <math.h>

using namespace cimg_library;
using namespace std;

//globa values
const double  epsilon=0.8; // the papamater smooth Dirac function (default value 1.5);
const float precision=0.009f;//precision of the error estimation

//functions
CImg<double> DiracU( CImg<double>& u0) ;
CImg<double> Heaviside(CImg<double>& u0);
CImg<double> ExtractContour(CImg<double> LevelSet);
CImg<unsigned char> get_level0(const CImg<>& img);

CImg<unsigned char> InitialLevelSet(CImg<double>&Img);
CImg<double> DiracF(CImg<double>& u1,CImg<double>& u2); 


//-----------------
// Main-MexFunction
//-----------------

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   if (nrhs < 13) mexErrMsgTxt("No enough input arguments.");
   if (nrhs >13) mexErrMsgTxt("Too many input arguments.");
   if (nrhs == 13){
         
       //Input Parameters ( inputs)
       CImg<double>  Img(prhs[0],true);         //input image
       CImg<double>  phi1(prhs[1],true);        //Rib cage curve
       CImg<double>  phi2(prhs[2],true);        //Second level set function which tracks the heart inside the rib cage
       CImg<int>     VolumeMask(prhs[3],true);
       const double dt= mxGetScalar( prhs[4]);     //time-step of iteration
       const double kappa = mxGetScalar(prhs[5]);  //coefficient of the weighted length term L(phi1)
       const double lambda1 = mxGetScalar(prhs[6]);//coefficient of insice C term
       const double lambda2 = mxGetScalar(prhs[7]);//coefficient of outside C term
       const double lambda = mxGetScalar(prhs[8]); //coefficient of the weighted length term L(phi2)
       const double mu = mxGetScalar(prhs[9]);     // coefficient of the internal (penalizing) energy term P(u)
       const double v = mxGetScalar(prhs[10]);     //coefficient of the weighted area term A(u)
       const unsigned int nb_iter = mxGetScalar(prhs[11]);//number of iterations
       CImg<double> g(prhs[12],true);              //edge indicator for phi2 evolution
       //////////////////////////////////////////////////////////////////////////////////////////////////
       
      //Design the initial distance functions phi1,phi2 
      unsigned char col[1]={2};//color filling
   
      phi1.draw_fill(0,0,col);
      phi2.draw_fill(0,0,col);
     
      //Define rib cage Mask;
      CImg<double> f(phi1.dimx(),phi1.dimy());
      f.fill(-2);
      cimg_forXY(phi1,x,y){if(phi1(x,y)==0) f(x,y)=1;}
     
      
      cimg_forXY(phi1,x,y){
        phi1(x,y)=-floor(phi1(x,y)-phi2(x,y))-1;
        phi2(x,y)=5*(phi2(x,y)-1);
      }
      phi1.distance_hamilton(15);//distance function
    //  CImgDisplay disp(phi2,"phi2",0);
        //Initializations
       CImg<double> dphi1(Img.dimx(),Img.dimy(),2); //derivatives of phi
       CImg<double> veloc1(phi1.dimx(),phi1.dimy());//evolution matrix
       CImg<double> N_dphi1(phi1.dimx(),phi1.dimy(),2); //Normalize gradient of function Phi1
       CImg<double> veloc2(phi2.dimx(),phi2.dimy());
       CImg<double> N_dphi2(phi2.dimx(),phi2.dimy(),2); //Normalize gradient of function Phi2
       
       //Chan-Vese Coefficients
       double c1=0, c2=0, Averagec1=0, Averagec2=0;
      
      //Edge indicator for Phi2 evolution
      CImg<double>dg(g.dimx(),g.dimy(),1,2);
      cimg_for3XY(g,x,y){
              dg(x,y,0)=0.5*(g(_n1x,y)-g(_p1x,y)), 
              dg(x,y,1)=0.5*(g(x,_n1y)-g(x,_p1y));
       }
              
      //Heaviside for the initial rib cage
     CImg<double> HeavisideF_R=Heaviside(f);
     
     double E1=1e20f,E2=1e20f;//Initial energies
     double Eold1 = 0, Eold2=0;  
     veloc1.fill(0);
     veloc2.fill(0);
     
     //////////////////////////////////////////////////////////////////////////////////////////////
     // PDEs
     for (unsigned int iter=0; iter<=nb_iter; iter++) {
                
                                       
             CImg<double>diracF1=DiracU(phi1);
             CImg<double>HeavisideF1=Heaviside(phi1);
             CImg<double>diracF2=DiracU(phi2);
             CImg<double>HeavisideF2=Heaviside(-phi2);
           
           //Estimation of derivatives of the Phi1,Phi2 and the chan-vese coefficients for the Phi1 evolution
            cimg_for3XY(phi1,x,y)if (VolumeMask(x,y)==0){
                
                 //Phi1-Chan Vese             
                 const double
                      phix1=0.5*(phi1(_n1x,y)-phi1(_p1x,y)),
                      phiy1=0.5*(phi1(x,_n1y)-phi1(x,_p1y));  //derivatives of phi1(central approximation)

                  const double Mag_dphi1= sqrt(pow(phix1,2)+ pow(phiy1,2)+1e-10); //magnitude of grad(phi)
                      N_dphi1(x,y,0)=phix1/Mag_dphi1;
                      N_dphi1(x,y,1)=phiy1/Mag_dphi1;

                 c1+=HeavisideF1(x,y)* Img(x,y);
                 c2+=(1-HeavisideF1(x,y))*Img(x,y);
                 Averagec1+=HeavisideF1(x,y);
                 Averagec2+=HeavisideF1(x,y);
                 /////////////////////////////////////////////////////////////////////////////////////
                 
                 //Phi2 evolution-Front propagation or Level set function without re-initialization
                 const double
                      phix2=0.5*(phi2(_n1x,y)-phi2(_p1x,y)),
                      phiy2=0.5*(phi2(x,_n1y)-phi2(x,_p1y));  //derivatives of phi2

                 const double Mag_dphi2= sqrt(pow(phix2,2)+ pow(phiy2,2)+1e-10); //magnitude of grad(phi2)
                      N_dphi2(x,y,0)=phix2/Mag_dphi2;
                      N_dphi2(x,y,1)=phiy2/Mag_dphi2;
             }

            //chan-vese coefficients update
            c1/=(Averagec1+1e-5);
            c2/=(Averagec2+1e-5);
                          
             
            cimg_for3XY(Img,x,y)if (VolumeMask(x,y)==0){
                     
                 //Chan-Vese level set function 
                const double
                    Laplac_phi1=(phi1(_n1x,y) + phi1(_p1x,y) + phi1(x,_n1y) + phi1(x,_p1y))-4*phi1(x,y),     //laplacian operator
                    K1=0.5*(N_dphi1(_n1x,y,0)-N_dphi1(_p1x,y,0))+0.5*(N_dphi1(x,_n1y,1)-N_dphi1(x,_p1y,1));//curvature estimation
                const double
                    phixx1=(phi1(_n1x,y)+phi1(_p1x,y)-2*phi1(x,y)),//second derivatives of Phi1
                    phiyy1=(phi1(x,_n1y)+phi1(x,_p1y)-2*phi1(x,y)); 
                
                //Evolution equation of Phi1
                veloc1(x,y)=mu*(Laplac_phi1-K1)-lambda1* diracF1(x,y)* pow(Img(x,y)-c1,2)+lambda2*diracF1(x,y)*pow(Img(x,y)-c2,2)+kappa*diracF1(x,y)*K1;
                E1+=lambda1*HeavisideF1(x,y)*pow(Img(x,y)-c1,2)+lambda2*(1-HeavisideF1(x,y))*pow(Img(x,y)-c2,2);  
                
                //Phi2-front propagation level set function (without re-initiallization)
                  const double
                     Laplac_phi2=(phi2(_n1x,y) + phi2(_p1x,y) + phi2(x,_n1y) + phi2(x,_p1y))-4*phi2(x,y),
                     K2=0.5*(N_dphi2(_n1x,y,0)-N_dphi2(_p1x,y,0))+0.5*(N_dphi2(x,_n1y,1)-N_dphi2(x,_p1y,1));
                  const double
                     phixx2=(phi2(_n1x,y)+phi2(_p1x,y)-2*phi2(x,y)),
                     phiyy2=(phi2(x,_n1y)+phi2(x,_p1y)-2*phi2(x,y));
               
                 
                  veloc2(x,y)=lambda* diracF2(x,y)*( dg(x,y,0)* N_dphi2(x,y,0) +dg(x,y,1)* N_dphi2(x,y,1) + g(x,y)*K2)+mu*(Laplac_phi2-K2)+v*g(x,y)*diracF2(x,y);
                  veloc2(x,y)=veloc2(x,y)*HeavisideF_R(x,y)*(1-HeavisideF1(x,y)*(-HeavisideF2(x,y)));
                //Energy estimation
                // E2+=lambda*g(x,y)*diracF2(x,y)* Mag_dphi2+1/2*mu*pow( Mag_dphi2-1,2)+v*HeavisideF2(x,y)*g(x,y);
//                   if (!(iter%400)) {
//                   get_level0(phi2).resize(disp.dimx(),disp.dimy()).draw_grid(20,20,0,0,false,false,col,0.4f,0xCCCCCCCC,0xCCCCCCCC).
//                   draw_text(5,5,"Iteration %d",col,0,1,11,iter).display(disp);
//                 }
               
            }
         
            phi1+=dt*veloc1;
            phi2+=dt*veloc2;
           
          if ((abs(Eold1-E1)<0.001f) && (abs(Eold2-E2)<0.001f)) break; 
          
            c1=0,Averagec1=0;
            c2=0,Averagec2=0;
            Eold1 = E1, Eold2=E2;
            E1=0;E2=0;
     
     }
  
    plhs[0]=  phi1.toMatlab();
    plhs[1]=  phi2.toMatlab();
  
    
  } 
    
   return;    
 
}
  

   
CImg<double> DiracU(CImg<double>& u0) {

  CImg<double> u(u0.dimx(),u0.dimy());
  u.fill(0);
   
  cimg_forXY(u0,x,y) { 
       if (u0(x,y)<=epsilon && u0(x,y)>=-epsilon){
           u(x,y)=(double)1/(2*epsilon)*(1+cos(3.14*u0(x,y)/epsilon));
       }
  
   }
   return u;
}

CImg<double> Heaviside(CImg<double>& u0) {

  CImg<double> u(u0.dimx(),u0.dimy());
  u.fill(0);
/*cimg_forXY(u0,x,y){
    u(x,y)=1/2*(1+2/3.14*atan(u(x,y)/epsilon));
}*/  
 cimg_forXY(u0,x,y) { 
                  
       if (u0(x,y)>=-epsilon && u0(x,y)<=epsilon){
           u(x,y)=(double) 1/2+u0(x,y)/(2*epsilon)+1/(2*3.14)*sin(3.14*u0(x,y)/epsilon);
       }
       if (u0(x,y)>epsilon) u(x,y)=1;
  
   }
    return u;
} 


/*******************************************************************************/
CImg<double> ExtractContour(CImg<double> LevelSet)
{
 CImg<double> Contour(LevelSet.dimx(),LevelSet.dimy(),1,1);
 Contour.fill(0);

 CImg_3x3(I,double);
 cimg_for3x3(LevelSet,x,y,0,0,I)
 {
  if(Icc*Icp<=0 || Icc*Icn<=0 || Icc*Ipc<=0 || Icc*Inc<=0)
   Contour(x,y) = 1;
 }
 return Contour;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Create a user-defined closed curve (Initial level set fuction)
CImg<unsigned char> InitialLevelSet(CImg<double>&Img){
       CImg<unsigned char> curve(Img.dimx(),Img.dimy(),Img.dimz(),2,0);
       unsigned char col1[2]={0,255}, col2[2]={200,255}, col3[2]={255,255};//colors
       curve.draw_grid(20,20,0,0,false,false,col1,0.4f,0xCCCCCCCC,0xCCCCCCCC).
       draw_text(5,5,"Please draw your curve\nin the middle of this window\n(Use your mouse)\n-heart initial curve",col1);
      CImgDisplay disp(curve,"Image",0);
       CImg<double> tempImg(Img);

       int xo=-1,yo=-1,x0=-1,y0=-1,x1=-1,y1=-1;
       while (!disp.is_closed && (x0<0 || disp.button)) {
        if (disp.button && disp.mouse_x>=0 && disp.mouse_y>=0) {
             if (x0<0) { xo = x0 = disp.mouse_x; yo = y0 = disp.mouse_y; } else {
                 x1 = disp.mouse_x; y1 = disp.mouse_y;
                 curve.draw_line(x0,y0,x1,y1,col2);//.display(disp);
                
                  tempImg.draw_point(x1,y1,col1).display(disp);
                 x0 = x1; y0 = y1;
              }
         }
         disp.wait();
        if (disp.is_resized) disp.resize(disp);
       }
 curve.draw_line(x1,y1,xo,yo,col2).channel(0).draw_fill(0,0,col3);
return curve;
}
//////////////////////////////////////////////////////////////////////////////////////////////




// get_level0() : Retrieve the curve corresponding to the zero level set of the distance function
//-------------
CImg<unsigned char> get_level0(const CImg<>& img) {
  CImg<unsigned char> dest(img);
  CImg_2x2(I,float); Inn = 0;
  cimg_for2x2(img,x,y,0,0,I) if (Icc*Inc<0 || Icc*Icn<0) dest(x,y) = 255; else dest(x,y) = Icc<0?100:0;
  return dest;
}
        
   
  
     

  


   
