function k_init = get_k_init(sensor_file_name)

euler = load(sensor_file_name);
euler(:,1) = [];

% focal length
f = 4.26;

points = f*[euler(:,2) euler(:,3)];
points(:,1) = points(:,1) - points(1,1);
points(:,2) = points(:,2) - points(1,2);

points2 = round(points);
xq=points2(:,1);
vq=points2(:,2);

% move to the center
ksize = 25;
% ksize = length(xq);
offset = [round(mean(xq)-ksize/2),round(mean(vq)-ksize/2)];
xq = xq - offset(1);
vq = vq - offset(2);

k = zeros(ksize);

for i=1:length(xq)
    k(vq(i),xq(i)) = 1;
end
k = k./sum(k(:));

k_init = rot90(k,-1);

% figure,imagesc(k)
% figure,imagesc(rot90(k,2))

end