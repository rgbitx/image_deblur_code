clear
close all

alldata = load('alldata.txt');
alldata(:,1)=[];

euler = load('sensorData.txt');
euler(:,1) = [];

% focal length
f = 4.26;

points = f*[euler(:,2) euler(:,3)];
points(:,1) = points(:,1) - points(1,1);
points(:,2) = points(:,2) - points(1,2);
points2 = 10*points;
pointsX = points2(:,1);
pointsY = points2(:,2);
rg_pointsX = range(pointsX);
rg_pointsY = range(pointsY);

% interp
x=[];
v=[];
flag = 1;
if rg_pointsX > rg_pointsY && rg_pointsX > 1
    x = pointsX;
    v = pointsY;
    flag = 1;
else
    x = pointsY;
    v = pointsX;
    flag = 2;
end
xq = ceil(min(x)):round(max(x));
vq0 = interp1(x,v,xq);
vq = round(interp1(x,v,xq));

% move to the center
ksize = 100;
% ksize = length(xq);
offset = [round(mean(xq)-ksize/2),round(mean(vq)-ksize/2)];
xq = xq - offset(1);
vq = vq - offset(2);

k = zeros(ksize);
if flag == 1
    for i=1:length(xq)
        k(ksize-xq(i)+1,vq(i)) = 1;
    end
    k = k./sum(k(:));
end
if flag == 2
   for i=1:length(xq)
        k(ksize-vq(i)+1,xq(i)) = 1;
   end 
   k = k./sum(k(:));
end

k1 = imresize(k,[25,25]);

figure;plot(xq,vq,'o')
figure;imagesc(k);
figure;imagesc(k1);

PSF = k1;
estimated_nsr = 0;
blurred = imread('image9.png');
wnr = deconvwnr(blurred, PSF, estimated_nsr);
figure; imshow(wnr);

% kernelSize = 25;
% points2 = points2*24/range(points2(:));
% minPoint = min(points2(:));
% points2 = points2 - minPoint;
% 
% k = zeros(kernelSize);
% 
% % set to center
% points3 = points2;
% xrang = range(points2(:,1));
% yrang = range(points2(:,2));
% 
% offset = [round(xrang/2-kernelSize/2),round(yrang/2-kernelSize/2)];
% points3(:,1) = points3(:,1) + offset(1);
% points3(:,2) = points3(:,2) + offset(2);
% 
% for i=1:length(points(:,1))
%    k(round(points3(i,1))+2,round(points3(i,2))+2) = 1;
% end
% 
% k = k/sum(k(:));

% figure
% imshow(k)


