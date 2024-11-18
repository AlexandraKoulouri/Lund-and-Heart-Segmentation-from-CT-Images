function [x,y] = MaxLikelihood(Img,RibCage,Displacements);

%Maximum likelihood to belong in the inner soft tissue
%The features are:
%Pixels Intensity
%Displacement values towards x and y axis


%define the region inside the rib cage
Mask=imfill(RibCage);
Mask=imerode(Mask, strel(ones(10,10)));

%find the pixels with gray intensities inside the rib cage
Pixels= find(Mask==1 & Img>0);

%Find the intensity values of the non-zero pixels
IntensityValues=Img(Pixels);

%Find the displacements values of the gray pixels
DisplacementXY=Displacements([Pixels,Pixels+size(Displacements,1)*size(Displacements,2)]);
DisplacementXY=reshape(DisplacementXY,[ length(Pixels),2]);


[Pmax,ClassNum,Pixels]= NaiveBayesianClassifier(Pixels,IntensityValues, DisplacementXY);
%ClassNum==1 means that the pixel belongs to heart, otherwise to
%background.

ProbabilitiesHeartMap=zeros(size(Img));
ProbabilitiesHeartMap(Pixels(ClassNum==1))=Pmax(ClassNum==1);
ProbabilitiesHeartMap = imerode(ProbabilitiesHeartMap ,strel(ones(10,10)));
Pixel = find(ProbabilitiesHeartMap==max(max(ProbabilitiesHeartMap)));

x=0;y=0;
%coordinates of the pixel
[y,x] = ind2sub([size(Img,1),size(Img,2)],Pixel(1));



