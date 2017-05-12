function [confidence,offset,scale] = calconfidence(gradient, ksize, step)

    [rows,cols] = size(gradient{1}); 
    
    range = floor(ksize/2);
    
    width = ksize;
    
    % pixAxis = offset + scale*(confAxis-1)
    offset = width;
    scale = step;
    % step of moving windows
%     step = floor(ksize/3);
    
    % save all confidence values
    % confidence row
    cfrow = 1;
    cfrows = ceil((rows-2*width)/step);
    cfcols = ceil((cols-2*width)/step);
    
    confidence = zeros(cfrows,cfcols);
    
    for row = width:step:(rows-width)
       cfcol = 1;
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
                        
           confidence(cfrow,cfcol)= sqrt(sumx^2+sumy^2)/(innerSum + 0.5);
           cfcol = cfcol + 1;
       end
       cfrow = cfrow + 1;
    end
    
end