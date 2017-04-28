function edge_weight = determine_truck(I_X, I_Y, I_mag,k1,k2)

% parameters
para = zeros(1,4);
num_angles = 4;
alpha = 9; %% parameter for attenuation of angles (must be odd)

angle_step = 2 * pi / num_angles;
angles = 0:angle_step:2*pi;
angles(num_angles+1) = []; % bin centers

I_theta = atan2(I_Y,I_X);%gradient direction
I_theta(find(isnan(I_theta))) = 0; % necessary????
[hgt, wid] = size(I_X);
% make orientation images
I_orientation = zeros([hgt, wid, num_angles]);

% for each histogram angle
cosI = cos(I_theta);
sinI = sin(I_theta);
for a=1:num_angles
    % compute each orientation channel
    tmp = (cosI*cos(angles(a))+sinI*sin(angles(a))).^alpha;
    tmp = tmp .* (tmp > 0);

    % weight by magnitude
    I_orientation(:,:,a) = tmp .* I_mag;
    t =  I_orientation(:,:,a);
    tt= sort(t(:),'descend');
    para(a) = tt(2*min(k1,k2));
    %para(a) = tt(2*sqrt(k1*k2));
end
edge_weight =min(para);
