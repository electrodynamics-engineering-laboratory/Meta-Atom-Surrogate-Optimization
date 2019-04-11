function Data= datainput_convex_I
%convex problem
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
Data.xlow=-10*ones(1,8); %variable lower bounds
Data.xup=10*ones(1,8); %variable upper bounds
Data.dim=8; %problem dimension
Data.integer=(1:8); %indices of integer variables
Data.continuous = []; %indices of continuous variables

%define objective function
Data.objfunction=@(x) 3.1*x(:,1).^2 + 7.6* x(:,2).^2 +6.9*x(:,3).^2 +0.004*x(:,4).^2 +...
    +19*x(:,5).^2 +3*x(:,6).^2 +x(:,7).^2  +4*x(:,8).^2 ;
end %function