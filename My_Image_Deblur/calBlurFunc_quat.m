
quat = load('sensorData.txt');
quat(:,1) = [];
rotMat = quat2rotm(quat);
[yaw,pitch,roll] = quat2angle(quat);
angles=[yaw,pitch,roll];
% focal length
f = 4.26;
K=[f 0 0
   0 f 0
   0 0 1];
K_inv = inv(K);

x0=0; y0=0;
point0 = [x0 y0 1]';

[row,col,num] = size(rotMat);
points = zeros(num,3);

for i=1:num
   H = K*rotMat(:,:,i)*K_inv;
   points(i,:) = H*point0;
end

points(:,1) = points(:,1)-points(1,1);
points(:,2) = points(:,2)-points(1,2);
points(:,3) = points(:,3)-points(1,3);

figure
plot(points(:,1),points(:,2))


