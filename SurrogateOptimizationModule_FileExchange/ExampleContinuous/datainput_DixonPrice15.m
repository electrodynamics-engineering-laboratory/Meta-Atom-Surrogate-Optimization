function Data = datainput_DixonPrice15

%15-dimensional dixon-price function
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

Data.xlow=-15*ones(1,15); %variable lower bounds
Data.xup=30*ones(1,15); %variable upper bounds
Data.dim=15; %problem dimension
%objective function
Data.objfunction=@(x) (x(1)-1)^2 + sum((2:Data.dim).*(2*x(2:Data.dim).^2-x(1:Data.dim-1)).^2,2);
Data.integer = []; %indices of integer variables
Data.continuous=(1:15); %indices of continuous variables
end %function