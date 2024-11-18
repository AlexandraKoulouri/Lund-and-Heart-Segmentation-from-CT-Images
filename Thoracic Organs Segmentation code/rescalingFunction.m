
function image = rescalingFunction(image,i_lower,i_upper);

image=double(image);
maxValue=double(i_upper);
i_upper=0.9*double(i_upper);i_lower=double(i_lower);
image(image>i_upper) = maxValue;
image(image<i_lower) = 0;
indeces = find(image>i_lower & image<i_upper);
image(indeces) = (double(image(indeces))-i_lower)/(i_upper-i_lower)*maxValue;
%image = medfilt2 ( image, [5,5] );
image=uint16(image);