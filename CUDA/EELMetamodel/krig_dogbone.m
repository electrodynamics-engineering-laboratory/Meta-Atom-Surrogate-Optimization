% Load data
raw_data = csvread('Dogbone.csv',1,0);

pl = raw_data(:,2);     % element length
pw = raw_data(:,1);     % element width
al = raw_data(:,3);     % arm length
aw = raw_data(:,4);     % arm width
cl = raw_data(:,5);     % core length
cw = raw_data(:,6);     % core width

mc = raw_data(:,7);     % magnitude co-pol
pc = raw_data(:,8);     % phase co-pol
mx = raw_data(:,9);     % magnitude x-pol
px = raw_data(:,10);    % phase x-pol

N = length(pl);         % number of data points

theta = 0.3*6;
a = 0;
sig2 = 1;

p = 2;

% Normalize data
d_max = 20;
d_min = 1;

pl = (pl - d_min)./(d_max - d_min);
pw = (pw - d_min)./(d_max - d_min);
al = (al - d_min)./(d_max - d_min);
aw = (aw - d_min)./(d_max - d_min);
cl = (cl - d_min)./(d_max - d_min);
cw = (cw - d_min)./(d_max - d_min);

% Pick modeling data
random_pick = 1;
pick_ratio = 0.25;

Np = round(pick_ratio * N);

if random_pick ~= 0
    pick_idx = randsample(N,Np);
else
    pick_idx = 1:Np;
end

% Pick test sites
test_idx = 1:N;
test_idx(pick_idx) = [];

%%

p=6;

% Compute distances
Sij = zeros(Np);
for ix = 1:Np
    for iy = 1:Np
        Sij(ix,iy) = (pl(pick_idx(ix)) - pl(pick_idx(iy)))^2 + (pw(pick_idx(ix)) - pw(pick_idx(iy)))^2 + ...
            (al(pick_idx(ix)) - al(pick_idx(iy)))^2 + (aw(pick_idx(ix)) - aw(pick_idx(iy)))^2 + ...
            (cl(pick_idx(ix)) - cl(pick_idx(iy)))^2 + (cw(pick_idx(ix)) - cw(pick_idx(iy)))^2;
    end
end

S0j = zeros(Np,N-Np);
for ix = 1:Np
    for iz = 1:(N-Np)
        S0j(ix,iz) = (pl(pick_idx(ix)) - pl(test_idx(iz)))^2 + (pw(pick_idx(ix)) - pw(test_idx(iz)))^2 + ...
            (al(pick_idx(ix)) - al(test_idx(iz)))^2 + (aw(pick_idx(ix)) - aw(test_idx(iz)))^2 + ...
            (cl(pick_idx(ix)) - cl(test_idx(iz)))^2 + (cw(pick_idx(ix)) - cw(test_idx(iz)))^2;
    end
end


% Compute covariances
Cij = (sig2 - a)*exp(-theta.*Sij);
C = padarray(Cij,[1 1],1,'post');
C(end,end) = 0;
C_inv = inv(C);

C0j = (sig2 - a)*exp(-theta*S0j);
c0 = padarray(C0j,[1 0],1,'post');

% Compute weights
w = zeros(Np+1,N-Np);
for iz =1:N-Np
    w(:,iz) = C_inv*c0(:,iz);
end

% Estimate magnitude
mc_est = zeros(N-Np,1);
mc_est_id = zeros(N-Np,1);
for ix = 1:N-Np
    mc_est(ix) = w(1:Np,ix)' * mc(pick_idx);
    if(find(~S0j(:,ix)))
        [ac dc] = find(~S0j(:,ix),1,'first');
        mc_est_id(ix) = mc(pick_idx(ac));
    else
        mc_est_id(ix) = (mc(pick_idx)'*(1./S0j(:,ix).^p))./sum(1./S0j(:,ix).^p);
    end
end

% Estimate phase
pc_est = zeros(N-Np,1);
pc_est_id = zeros(N-Np,1);
for ix = 1:N-Np
    pc_est(ix) = w(1:Np,ix)' * pc(pick_idx);
    if(find(~S0j(:,ix)))
        [ac dc] = find(~S0j(:,ix),1,'first');
        pc_est_id(ix) = pc(pick_idx(ac));
    else
        pc_est_id(ix) = (pc(pick_idx)'*(1./S0j(:,ix).^p))./sum(1./S0j(:,ix).^p);
    end
end

% Estimate root mean errors
mc_er = sqrt(sum((mc(test_idx) - mc_est).^2)./(N-Np));  % Magnitude error
mc_id_er = sqrt(sum((mc(test_idx) - mc_est_id).^2)./(N-Np));
pc_er = sqrt(sum((pc(test_idx) - pc_est).^2)./(N-Np));  % Phase error
pc_id_er = sqrt(sum((pc(test_idx) - pc_est_id).^2)./(N-Np));

[Np mc_er pc_er mc_id_er pc_id_er]
[mc_id_er pc_id_er];

