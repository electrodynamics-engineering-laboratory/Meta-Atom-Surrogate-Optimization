function Data= datainput_GWSS20_I
%--------------------------------------------------------------------------
%Copyright (c) 2012 by Juliane Mueller
%
% This file is part of the surrogate model module toolbox.
%
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%Tampere University of Technology, Finland
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
Data.xlow=zeros(1,20); %variable lower bounds
Data.xup=ones(1,20); %variable upper bounds
%define data needed for computing objective function
StatPoint = [...
0.1113  0.4312  0.4005  0.0290  0.3407  0.7404  0.3254  0.8732  0.2982  0.4647  0.5835  0.9671  0.6511  0.7339  0.7554  0.7001  0.2046  0.4895  0.0168  0.2596  
0.0676  0.2388  0.3520  0.2962  0.5106  0.8120  0.9443  0.4184  0.1276  0.2262  0.9370  0.4063  0.8412  0.0365  0.8183  0.6819  0.7485  0.2394  0.8701  0.5953  
0.7606  0.7752  0.2521  0.9394  0.8989  0.4783  0.3207  0.6626  0.7606  0.1632  0.4409  0.8366  0.4924  0.4403  0.9590  0.0166  0.3428  0.7571  0.7471  0.6684  
0.3112  0.0075  0.1047  0.0416  0.0151  0.7655  0.7663  0.6837  0.1840  0.0986  0.7070  0.8130  0.9253  0.0989  0.0166  0.1438  0.1840  0.2523  0.4117  0.6720  
0.8173  0.5537  0.4085  0.2906  0.2303  0.4465  0.6033  0.0256  0.3248  0.4831  0.9061  0.7732  0.2874  0.7148  0.4714  0.3953  0.5526  0.4326  0.5315  0.8252  
0.1549  0.2473  0.5034  0.3881  0.6674  0.7438  0.1361  0.6293  0.5676  0.4247  0.4787  0.5082  0.4753  0.4324  0.0017  0.7740  0.9442  0.4576  0.8027  0.4941];  
StatValue = [...
0.0000
9.9836
0.0000
-21.1425
44.1928
-49.1906];
%objective function handle
Data.objfunction=@(x)GWSS20_2_11(x,StatPoint,StatValue);
Data.dim=20; %problem dimension
Data.integer =(1:20); %indices of integer variables
Data.continuous = []; %indices of continuous variables
end %function

function y=GWSS20_2_11(x,StatPoint,StatValue) %objective function 
    k = length(StatValue);
    y=zeros(size(x,1),1);
    for ii = 1:size(x,1)
        Xtemp = repmat(x(ii,:),k,1);
        TempNorm = sum((Xtemp-StatPoint).^2,2);
        TempProd = ones(1,k);
        for i = 1:k
            TempProd(i) = prod(TempNorm(1:(i-1)))*prod(TempNorm((i+1):k));
        end
        y(ii,1) = (TempProd*StatValue)/sum(TempProd);
    end
end %GWSS20_2_11 