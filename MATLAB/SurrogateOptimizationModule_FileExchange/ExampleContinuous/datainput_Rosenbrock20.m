function Data = datainput_Rosenbrock20
%20-dimesnional Rosenbrock function
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

Data.xlow=-5*ones(1,20); %variable lower bounds
Data.xup=10*ones(1,20); %variable upper bounds
Data.dim = 20; %problem dimension
%objective function
Data.objfunction= @(x) sum(100*(x(1:Data.dim-1).^2-x(2:Data.dim)).^2 + (x(1:Data.dim-1)-1).^2,2);
Data.integer=[]; %indices of integer variables
Data.continuous = (1:20);  %indices of continuous variables

end %function