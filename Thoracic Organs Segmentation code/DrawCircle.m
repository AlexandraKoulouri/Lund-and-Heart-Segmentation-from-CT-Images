function [Img] = DrawCircle(R,c1,c2,x,y,lineWidth);
Img = zeros(x,y);
theta = [0 : 0.01 : 2 * pi];
stack1 = 0;
if nargin == 5
    
    lineWidth = 1;
end
for i = 1 : length( theta )

        temp_x = 0; temp_y = 0;   
             
        temp_x =ceil( R* cos (theta(i)) + c1);
        
        temp_y = ceil (R* sin (theta(i)) + c2);
        stack1 = stack1 + 1;
        I_new(stack1) = temp_x ;
        J_new(stack1) = temp_y  ;
        Img( temp_x , temp_y ) = 1;
       
        %plot(I_new,J_new,'.',c1,c2,'r*')
end

se1 = strel('line',lineWidth,0);
se2 = strel('line',lineWidth,90);
Img = imdilate(Img,[se1 se2],'same');
