function Data= datainput_ex_I
%from http://www.aridolan.com/ga/gaa/MultiVarMin.html
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
%define objective function
Data.objfunction=@(x) x(:,1).*sin(x(:,1)) + 1.7*x(:,2).*sin(x(:,1)) - ...
    1.5*x(:,3) - 0.1*x(:,4).*cos(x(:,4)+x(:,5)-x(:,1)) +...
    (0.2*x(:,5).^2-x(:,2)) - 1;
%no constraints
Data.xlow=-100*ones(1,5);%variabe lower bounds
Data.xup=100*ones(1,5); %variabe upper bounds
Data.dim= 5; %problem dimension
Data.integer=(1:5); %indices of integer variables
Data.continuous = []; %indices of continuous variables
end %function