function [f,g] = bdc_obj(k,X,flipX,Y,mode)

switch mode
   case 'psf'
      [f,g] = kernel_obj(k,X,flipX,Y);
   case 'img'
      [f,g] = image_obj(k,X,flipX,Y);
end

end

function [f,g] = image_obj(k, X, flipX, Y)
  
res = {conv2(X{1}, k, 'same') - Y{1}, ...
       conv2(X{2}, k, 'same') - Y{2}};
    
f = 1/2*(norm(vec(res{1}))^2 + norm(vec(res{2}))^2);
g = [vec(conv2(res{1},rot90(k,2),'same')); vec(conv2(res{2},rot90(k,2),'same'))];

end


function [f,g] = kernel_obj(k, X, flipX, Y)
  
f = 0;
g = zeros(size(k));
khs = floor(size(k,1)/2);

for i = 1:length(X)
  r = conv2(X{i}, k, 'same')-Y{i};
  f = f+1/2*norm(r,'fro')^2;
  g = g + conv2(cpad(r,khs),flipX{i},'valid');
end  
   
g = vec(g);
end


function y = cpad(x,khs)
   y = zeros(size(x)+2*khs);
   y(khs+1:end-khs,khs+1:end-khs) = x;
end
