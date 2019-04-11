function Data = datainput_Shubert
%2-dimensional Shubert function
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

Data.xlow = [-10, -10]; %lower variable bounds
Data.xup = [10,10]; %upper variable bounds
Data.dim=2; %problem dimesnion
Data.integer=[]; %indices of integer variables
Data.continuous=[1,2]; %indices of continuous variables
I1=(1:5)';
I2=I1+1;
%objective function
Data.objfunction = @(x) sum((I1.* cos(I2.*x(1)+I1)))*sum((I1.* cos(I2.*x(2)+I1)));

end %function