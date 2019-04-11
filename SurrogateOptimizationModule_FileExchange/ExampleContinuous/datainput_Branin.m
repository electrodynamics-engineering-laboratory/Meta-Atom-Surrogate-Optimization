function Data = datainput_Branin
%Branin function
%3 global minima, no local  minima
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

Data.xlow=[-5,0]; %variable lower bounds
Data.xup=[10,15]; %variable upper bounds
%objective function
Data.objfunction= @(x)(x(:,2)-5.1*x(:,1).^2./(4*pi^2)+5*x(:,1)./pi-6).^2 + 10*(1-1/(8*pi))*cos(x(:,1))+10;
Data.integer=[]; %indices of integer variables
Data.continuous=(1:2); %indices of continuous variables
Data.dim = 2; %problem dimension
end %function