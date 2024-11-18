function [Thres] = Histogr(Img);

Totalnum = size(Img,1)*size(Img,2)*size(Img,3);
Img=reshape(Img, [Totalnum,1]);
[Num,binsValue]=hist(double(Img),100);

DesireNum = round(0.9*Totalnum);

Cum=zeros(length(Num),1);
Cum(1)=Num(1);

for i=2:length(Num)
    Cum(i)=Cum(i-1)+Num(i);
end

i=find(Cum>=0.955*max(Cum));
plot(binsValue,Cum)
Thres=binsValue(i(1));

