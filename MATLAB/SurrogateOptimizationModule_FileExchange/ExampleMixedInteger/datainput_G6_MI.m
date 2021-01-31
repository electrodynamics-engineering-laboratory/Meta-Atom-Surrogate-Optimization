function Data=datainput_G6_MI

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
%objective function handle
Data.objfunction=@(x) (x(:,1)-10).^3+(x(:,2)-20).^3;
%constraint function handles
Data.constraint{1}=@(x) -(x(:,1)-5).^2 -(x(:,2)-5).^2+100;
Data.constraint{2}=@(x) (x(:,1)-6).^2+(x(:,2)-5).^2-82.81;
Data.xlow=[13,0]; %variable lower bounds
Data.xup=[100,100]; %variable upper bounds
Data.integer=1; %indices of integer variables
Data.continuous=2; %indices of continuous variables
Data.dim = 2; %problem dimension
end %function