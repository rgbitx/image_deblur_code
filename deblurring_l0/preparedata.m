function [B_org,B] = preparedata(img_name)

    B_org = imread(img_name);
    
    if size(B_org,3)==3   
        B = rgb2gray(B_org);
        B = double(B)/255;
    else
        B = double(B_org)/255;        
    end
end