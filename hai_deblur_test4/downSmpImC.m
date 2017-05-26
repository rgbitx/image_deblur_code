function sI=downSmpImC(I,ret)
% refer to Levin's code
if (ret==1)
    sI=I;
    return
end


sig=1/pi*ret;

g0=[-50:50]*2*pi;
sf=exp(-0.5*g0.^2*sig^2);
sf=sf/sum(sf);
csf=cumsum(sf);
csf=min(csf,csf(end:-1:1));
ii=find(csf>0.05);

sf=sf(ii);
sum(sf);

I=conv2(sf,sf',I,'valid');

[gx,gy]=meshgrid([1:1/ret:size(I,2)],[1:1/ret:size(I,1)]);

sI=interp2(I,gx,gy,'bilinear');
