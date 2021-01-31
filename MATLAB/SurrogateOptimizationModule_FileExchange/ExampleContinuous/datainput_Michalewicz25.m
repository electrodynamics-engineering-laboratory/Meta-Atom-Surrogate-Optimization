function Data = datainput_Michalewicz25
%25-dimensional Michalewicz function
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

Data.xlow=zeros(1,25);  %lower variable bounds
Data.xup=pi*ones(1,25); %upper variable bounds
Data.dim = 25; %problem dimension
Data.objfunction= @(x) -sum(sin(x).*(sin((1:Data.dim).*x.^2/pi)).^20,2); %objective function
Data.integer = []; %indices of integer variables
Data.continuous = (1:25); %indices of continuous variables


end %function