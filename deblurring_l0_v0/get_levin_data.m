function B = get_levin_data()

load('im05_flit01.mat');
K1_gt = f;% ground true kernel
L_gt = x;% ground true image
B1 = y;
load('im05_flit03.mat');
K2_gt = f;% ground true kernel
B2 = y;
%% input
B = zeros([size(B1),2]);% blurry image
B(:,:,1) = B1;
B(:,:,2) = B2;   

end