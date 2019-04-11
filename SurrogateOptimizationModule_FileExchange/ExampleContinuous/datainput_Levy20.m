function Data = datainput_Levy20

%20-dimensional Levy function
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

Data.xlow=-10*ones(1,20); %lower variable bounds
Data.xup=10*ones(1,20); %upper variable bounds
Data.dim = 20; %problem dimesnion
y=@(x) 1+(x-1)/4;
Data.objfunction=@(y) sin(pi*y(1))^2 + sum((y(1:Data.dim-1).^2.*(1+10*sin(pi*y(1:Data.dim-1)).^2)),2)...
    +(y(Data.dim)-1)^2*(1+sin(2*pi*y(Data.dim))^2);
Data.integer =[]; %indices of integer variables
Data.continuous = (1:20); %indices of continuous variables

end %function