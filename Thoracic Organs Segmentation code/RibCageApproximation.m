function [FinalCurve]=RibCageApproximation(I,Thres);
%Draw a  closed curve which approximates the cross
%section of the rib cage in CT slices.
%      [FinalCurve]=RibCageApproximation(I,Thres);
%       Input parameters
%       Img=> the CT slice
%       Thres=> the desired threshold for the ribs (bones) extraction.

Img=squeeze(I);

%apply  threshold
BinaryImage=Img>Thres;
%perform morphological open operation
BinaryImage = imopen(BinaryImage,strel([1 1;1 1]));

%find the pixels of the rib cage.
[x,y]=find(BinaryImage==1);

%define the center of the image.
ImgCenter =round( [size(BinaryImage)]./2);


%Estimate the vector between bone pixels and image centre
vector = [x-ImgCenter(1),y-ImgCenter(2)];
%Vector's angle in degrees.
Angle = mod(360+180/pi*atan2(vector(:,2),vector(:,1)),360);


%In each directions (0, 10,...., 360 degrees), define the bone pixel which is closer to the image
%center and keep it for the interpolation process.

Theta=[0:10:360];k=0;

Ind=cell([length(Theta)-1,1]);
for i =1:length(Theta)-1
    [Ind{i,1}] = find(Angle>=Theta(i) & Angle<Theta(i+1));
    if ~isempty(Ind{i,1})
      IndMinDist = FunctionMinDistance(x(Ind{i,1}),y(Ind{i,1}),ImgCenter(1),ImgCenter(2));
      k=k+1;
      IndOfMinDistance(k,1)=Ind{i,1}(IndMinDist(1),1);
       
    end 
end


%Set points with the minimum distance from the center of the ribs cavity in an order
[xi, yi]=PixelsInOrder(x(IndOfMinDistance),y(IndOfMinDistance),ImgCenter(1),ImgCenter(2));

%Cubic interpolation
[xii, yii] = Curvefitting(xi,yi,0.2);
xii = round(xii);yii=round(yii);

%Draw the final curve
FinalCurve=zeros(size(BinaryImage));

for i =1:length(xii)
    FinalCurve(xii(i),yii(i))=1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function  IndOfMinDistance = FunctionMinDistance(I,J,c1,c2)
%This function calculates the distance between each contour (I,J) point and
%the centre of the contour
%and finds the edge with the  minimum distance 
%      IndOfMinDistance = FunctionMinDistance(I,J,c1,c2)

%      Input Parameters
%      [I,J] points coordinates
%      [c1,c2] coordinates of the centre of the contour

MinDistance = zeros(length(I),1);

MinDistance(:) = sqrt((I(:)-c1).^2+(J(:)-c2).^2);
IndOfMinDistance =find (MinDistance <= min (MinDistance));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 function [Final_x, Final_y]=PixelsInOrder(x,y,c1,c2);

%Define a cartesian system xy with center[c1,c2],
%Set the points [x,y] in a circular order according to theirs angles 
%       [Final_x, Final_y]=PixelsInOrder(x,y,c1,c2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 c1=round(c1);c2=round(c2);

 vector = [x-c1,y-c2];
 %estimate the angle between the cartesian axis x and the points [x,y]
 Theta = mod(360+180/pi*atan2(vector(:,2),vector(:,1)),360);
 
 %sort the [x,y] points 
 [ThetaSorted,Ind] = sort(Theta);
 
 Final_x=zeros(length(x),1);
 Final_y=zeros(length(y),1);
 Final_x(1:length(x)) =x(Ind);
 Final_y(1:length(y)) =y(Ind);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [xi,yi] = Curvefitting(x,y,RES)
% Interpolate(cubic) the points of the image in order to create a close curve
%     [xi,yi] =Curvefitting(x,y,RES)
%
%      RES: resolution desired

%     
%      Modified Code by A. Koulouri
%      
%      Chenyang Xu and Jerry L. Prince, 3/8/96, 6/17/97
%      Copyright (c) 1996-97 by Chenyang Xu and Jerry L. Prince
%      Image Analysis and Communications Lab, Johns Hopkins University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert to column vector
x = x(:); y = y(:);

N = length(x);  

% make it a circular list since we are dealing with closed contour
 x = [x;x(1)];
 y = [y;y(1)];

dx = x([2:N+1])- x(1:N);
dy = y([2:N+1])- y(1:N);
d = sqrt(dx.*dx+dy.*dy);  % compute the distance from previous node for point 2:N+1

d = [0;d];   % point 1 to point 1 is 0 

% now compute the arc length of all the points  from  the first point [x(1),
% y(1)].
% we use matrix multiply to achieve summing 
M = length(d);
d = (d'*uppertri(M,M))';

% now ready to reparametrize the closed curve in terms of arc length
maxd = d(M);

if (maxd/RES<3)
   error('RES too big compare to the length of original curve');
end

di = 0:RES:maxd;

xi = interp1(d,x,di,'cubic');
yi = interp1(d,y,di,'cubic');

N = length(xi);

if (maxd - di(length(di)) <RES/2)  % deal with end boundary condition
   xi = xi(1:N-1);
   yi = yi(1:N-1);
end

function X = uppertri(M,N)
% UPPERTRI   Upper triagonal matrix 
%            UPPER(M,N) is a M-by-N triagonal matrix

[J,I] = meshgrid(1:M,1:N);
X = (J>=I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




