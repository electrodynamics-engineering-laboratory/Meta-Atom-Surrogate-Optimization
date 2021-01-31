function Data = datainput_Sphere27
%27-dimesnional Sphere (de Jong) function
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

Data.xlow=-5.12*ones(1,27);  %variable lower bounds
Data.xup=5.12*ones(1,27); %variable upper bounds
Data.dim = 27; %problem dimension
Data.integer=[];  %indices of integer variables
Data.continuous =(1:27); %indices of continuous variables
Data.objfunction= @(x) sum(x.^2,2); %objective function
end %function