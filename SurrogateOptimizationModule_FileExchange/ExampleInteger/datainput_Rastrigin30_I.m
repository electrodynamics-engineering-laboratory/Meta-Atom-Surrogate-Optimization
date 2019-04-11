function Data= datainput_Rastrigin30_I
%30-dimensional Rastrigin function
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
Data.xlow=-ones(1,30); %variable lower bounds
Data.xup=3*ones(1,30); %variable upper bounds
%objective function
Data.objfunction=@(x)sum((x.^2) - cos(2*pi*x),2);
%no constraints
Data.integer =(1:30); %indices of integer variables
Data.continuous=[]; %indices of continuous variables
Data.dim = 30; %problem dimension
end %function