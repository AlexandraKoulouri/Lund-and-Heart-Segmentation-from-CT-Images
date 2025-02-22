function [Pmax,ClassNum,dataPoints]= NaiveBayesianClassif(dataPoints,Intensity, Displacements);

%load mean and std values of features extracted by hand
load MeanValues.mat
load StdValues.mat
%the first column of the .mat files corresponds to the manually segmented heart and the
%second column to background values .

%Intensity 
Intensity=double(meshgrid(Intensity,1:2));
MeanValue=meshgrid(MeanValues(1,1:2),1:size(dataPoints,1))';
StdValue=meshgrid(StdValues(1,1:2),1:size(dataPoints,1))';
A1=exp(-(Intensity-MeanValue).^2./(2*StdValue.^2));
clear Intensity StdValue MeanValue


%Displacement x
DisplacementX=meshgrid(Displacements(:,1),1:2);
MeanValue=meshgrid(MeanValues(2,1:2),1:size(dataPoints,1));
StdValue=meshgrid(StdValues(2,1:2),1:size(dataPoints,1));
A2=exp(-(DisplacementX-MeanValue').^2./(2*StdValue'.^2));
clear DisplacementX StdValue MeanValue
%Displacement y
DisplacementY=meshgrid(Displacements(:,2),1:2);
MeanValue=meshgrid(MeanValues(3,1:2),1:size(dataPoints,1));
StdValue=meshgrid(StdValues(3,1:2),1:size(dataPoints,1));
A3=exp(-(DisplacementY-MeanValue').^2./(2*StdValue'.^2));
clear DisplacementY StdValue MeanValue


P=A1.*A2.*A3;
clear A1 A2 A3 

[Pmax,ClassNum]=max(P,[],1);
