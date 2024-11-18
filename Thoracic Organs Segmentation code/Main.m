%Before .m execution add paths
Path1= cd;  Path1=[Path1,'\OpticalFlow3D'];path(Path1,path);
Path1 = cd; Path1 = [Path1, '\MAP']; path(Path1, path);
Path1 = cd; Path1 = [Path1, '\TwoPhaseLevel Sets']; path(Path1, path);

%input images
%Inhale and Exhale  images
load InhaleImg.mat;load ExhaleImg.mat;

% InhaleImg=imresize(imread('InhaleImg.jpg'),0.3);
% ExhaleImg=imresize(imread('ExhaleImg.jpg'),0.3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Img=reshape([InhaleImg;ExhaleImg]',[size(ExhaleImg,2),size(InhaleImg,1),2]);

%Rib Cage Approximation
Thresh = Histogr(Img);


%MaxValue=max(max(max(Img)));
%Img = rescalingFunction(Img,0.16*MaxValue,MaxValue);
   
[RibCage]=RibCageApproximation(Img(:,:,1),Thresh);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Feature extraction 

%Optical Flow for the desplacement estimation
 i1= permute(squeeze(Img(:,:,1)),[1 2 3]);%Inhale Img
 i2= permute(squeeze(Img(:,:,2)),[1 2 3]);%Exhale Img
 Displacement = OpticalFlow3D(i1,i2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum likelikehood pixel to belong Heart.  
[J,I] = MaxLikelihood(Img(:,:,1),RibCage,squeeze(Displacement));
 
 

InnerCurve = DrawCircle(4,I(1),J(1),size(Img,1),size(Img,2),1);

%Level Sets Process
[phi1,phi2]=TwoPhasesLevelSet(Img(:,:,1),RibCage,InnerCurve);

figure
 imshow(Img(:,:,1),[]); hold on
 c1 = contour(phi1,[0 0],'r','LineWidth',2.2);
 c2 = contour(phi2,[0 0],'g','LineWidth',2);
 title('Zero Level Set')
 hold off;

  
 [I_heart]=find((phi2)<=0);
 [I_lung]=find((phi1)<=0);

 SegmentOfThorax= zeros(size(Img));
 SegmentOfThorax(I_heart)=1;
 SegmentOfThorax(I_lung)=1;
 SegmentOfThorax= imfill(SegmentOfThorax);

 figure(1)
 imshow(Img(:,:,1),[]); hold on
 c = contour(SegmentOfThorax,[0 0],'m','LineWidth',2);
 hold off

 
