function [nLevels,mbsz,mksz] = prepare_multi_scale(bsz,ksz,min_ksz,scale_factor)
%input: bsz: 2X1, for blurry image
%       ksz: 2X1, for kernel
% scale_factor:
%output:
% nLevels: n
% mbsz: nX2, for blurry image
% mksz: nX2, for kernel

if ~exist('min_ksz','var') || min_ksz < 3
    min_ksz = 3;
end
if ~exist('scale_factor','var')
    scale_factor = sqrt(2);
end

min_ksz = min_ksz+(1-mod(min_ksz,2));% odd number
rr = bsz(1); cc = bsz(2);
rr_k = ksz(1); rr_k = rr_k+(1-mod(rr_k,2));
cc_k = ksz(2); cc_k = cc_k+(1-mod(cc_k,2));
nLevels = ceil(-log(min_ksz/min(rr_k,cc_k))/log(scale_factor))+1;

mbsz = zeros(nLevels,2);
mksz = zeros(nLevels,2);
mksz(nLevels,:) = [rr_k,cc_k];
mbsz(nLevels,:) = [rr,cc];
iLevel = nLevels-1;
count = 1;
while min(mksz(iLevel+1,:)) > min_ksz
    mksz(iLevel,:) = round(mksz(iLevel+1,:)/scale_factor);    
    mksz(iLevel,:) = mksz(iLevel,:)+(1-mod(mksz(iLevel,:),2));%odd number
    if mksz(iLevel,:) == mksz(iLevel+1,:)
        mksz(iLevel,:) = mksz(iLevel,:)-2;
    end    
    if mksz(iLevel,1) < min_ksz
        mksz(iLevel,1) = min_ksz;        
    end
    if mksz(iLevel,2) < min_ksz
        mksz(iLevel,2) = min_ksz;       
    end
    % correct scale factor
    scale_factor_correct = mksz(iLevel+1,:)./mksz(iLevel,:);
    mbsz(iLevel,:) = round(mbsz(iLevel+1,:)./scale_factor_correct);
    mbsz(iLevel,:) = mbsz(iLevel,:)-(1-mod(mbsz(iLevel,:),2));%odd number 
    count = count+1;
    iLevel = iLevel-1;
end
nLevels = count;
mbsz = mbsz(end-count+1:end,:);
mksz = mksz(end-count+1:end,:);

end%function