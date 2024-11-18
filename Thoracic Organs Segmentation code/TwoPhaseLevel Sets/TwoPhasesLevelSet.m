function [phi1,phi2]=TwoPhasesLevelSet(Img,phi1,phi2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Img=double(Img);
Img = (double(Img/max(max(Img))*255));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%coefficients for Phi2 evolution
sigma=0.9;                               % scale parameter in Gaussian kernel for smoothing.
G=fspecial('gaussian',5,sigma);
Img_smooth=conv2(double(Img),G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2; 
g=1./(1+abs(f).^1);                      % edge indicator function.
dt= 10;
mu=0.19/dt;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define where is the background and where is the volume
VolumeMask=fillingProc(2,2,Img);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

[phi1,phi2]= EvolutionProcess(medfilt2(Img,[3,3]),phi1,phi2,VolumeMask,dt,16,.1,.1,6,mu,-3.5,1000,g);


