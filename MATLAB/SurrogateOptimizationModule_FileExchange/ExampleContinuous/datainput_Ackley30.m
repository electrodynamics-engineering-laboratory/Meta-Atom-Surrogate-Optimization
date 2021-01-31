function Data = datainput_Ackley30
%30-dim ackley function

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

Data.xlow=-15*ones(1,30); %lower variable bounds
Data.xup=30*ones(1,30); %upper variable bounds
Data.dim=30;    %problem dimension
%objective function
Data.objfunction=@(x)-20*exp(-0.2*sqrt(sum(x.^2,2)/Data.dim)) - exp(sum(cos(2*pi*x),2)/Data.dim);
Data.integer=[]; %indices of integer variables
Data.continuous = (1:30); %indices of continuous variables
end %function