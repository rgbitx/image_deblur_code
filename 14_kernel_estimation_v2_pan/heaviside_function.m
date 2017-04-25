function mask = heaviside_function(I, thresh,level)

          if (~exist('level','var'))
              level  =2;
          end
           I0 = zeros(size(I));
           I0(find(I<thresh)) = 0;
           I0(find(I>=thresh)) = 1;
           mask = I0;
%            [i,j] =find(I0==1);
%            %if length(i)>(size(I,1)*size(I,2)/2)
%            if length(i)>2*min(size(I,1),size(I,2))
%                mask = I0;
%            else
%                kk = sort(I(:),'descend');
%                thresh = kk(floor(length(kk)*0.5*thresh)+1);
%                I0(find(I<thresh)) = 0;
%                I0(find(I>=thresh)) = 1;
%                mask = I0;
%            end
end