function Data = datainput_Rosenbrock10
%10-dimensional Rosenbrock function

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

Data.xlow=-5*ones(1,10); %lower variable bounds
Data.xup=10*ones(1,10); %upper variable bounds
%objective function
Data.objfunction= @(x) sum(100*(x(1:9).^2-x(2:10)).^2 + (x(1:9)-1).^2,2);
Data.dim=10; %problem dimension
Data.integer =[]; %indices of integer variables
Data.continuous=(1:10); %indices of continuous variables
end %function