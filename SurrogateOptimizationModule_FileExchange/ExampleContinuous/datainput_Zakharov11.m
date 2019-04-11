function Data = datainput_Zakharov11
%11-dimensional Zakharov function
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

Data.xlow=-5*ones(1,11); %variable lower bounds
Data.xup=10*ones(1,11); %variable upper bounds
Data.dim = 11; %problem dimension
Data.integer=[]; %indices of integer variables
Data.continuous = (1:11); %indices of continuous variables
%objective function
Data.objfunction= @(x) sum(x.^2,2) + 1/4*sum(((1:Data.dim).*x).^2,2) + 1/16*sum(((1:Data.dim).*x).^4,2);
end %function