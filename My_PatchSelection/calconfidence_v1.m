function fp = calconfidence(gradient, ksize)

    [rows,cols] = size(gradient{1}); 
    
    range = floor(ksize/2);
    
    width = ksize;
    
    step = floor(width/3);
    
    fp = 0;
    
    cnt = 0;
    
    for row = width:step:(rows-width)
       for col = width:step:(cols-width)
           
           xmin = row-range;
           xmax = row+range;
           ymin = col-range;
           ymax = col+range;
           
           gradx = gradient{1}(xmin:xmax,ymin:ymax);
           grady = gradient{2}(xmin:xmax,ymin:ymax);
           
           sumx = sum(gradx(:));
           sumy = sum(grady(:));
           
           innerSqrt = sqrt(gradx.^2+grady.^2);
           innerSum = sum(innerSqrt(:));
           
           fp = fp + sqrt(sumx^2+sumy^2)/(innerSum + 0.5);
           cnt = cnt + 1;
       end
    end
    
    assert(cnt~=0,'fp count should be nonzero!');
    
    fp = fp/cnt;
    
end