function Data= datainput_Rastrigin30_MI

%unconstrained
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
Data.objfunction=@(x)sum((x.^2) - cos(2*pi*x),2); %objective function handle
Data.integer=(1:10); %indices of integer variables
Data.continuous=(11:30); %indices of continuous variables
Data.dim=30; %problem dimension
end %function