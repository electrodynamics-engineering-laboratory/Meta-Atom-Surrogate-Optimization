function Data= datainput_ex_MI
%from http://www.aridolan.com/ga/gaa/MultiVarMin.html
% no constraints

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
Data.xlow=-100*ones(1,5); %variable lower bounds
Data.xup=100*ones(1,5); %variable upper bounds
%objective function handle
Data.objfunction=@(x) x(:,1).*sin(x(:,1)) + 1.7*x(:,2).*sin(x(:,1)) - ...
    1.5*x(:,3) - 0.1*x(:,4).*cos(x(:,4)+x(:,5)-x(:,1)) +...
    (0.2*x(:,5).^2-x(:,2)) - 1;
Data.integer=(1:2); %indices of integer variables
Data.continuous=(3:5); %indices of continuous variables
Data.dim = 5;
end %function