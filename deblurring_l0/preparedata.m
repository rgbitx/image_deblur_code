function [B_org,B] = preparedata(img_name_set)
n = length(img_name_set);
B_org = cell(n,1);
for ii = 1:n
    B_org{ii} = imread(img_name_set{ii});
end
[rr,cc,dd] = size(B_org{1});
B = zeros(rr,cc,n,'uint8');
if dd > 1
    for ii = 1:n
        B(:,:,ii) = rgb2gray(B_org{ii});
    end
else
    for ii = 1:n
        B(:,:,ii) = B_org{ii};
    end    
end


end