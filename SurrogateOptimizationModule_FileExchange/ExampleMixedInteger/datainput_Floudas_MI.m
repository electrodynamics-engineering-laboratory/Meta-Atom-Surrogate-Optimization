function Data = datainput_Floudas_MI

%* MINLP literature problem
%* Floudas, 1995.

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
Data.xlow =[0,0.2,-2.22554]; %variable lower bounds
Data.xup=[1,1,-1]; %variable upper bounds
Data.dim=3;
Data.integer=1; %indices of integer variables
Data.continuous=(2:3); %indices of continuous variables

Data.objfunction = @(x) 5*(x(:,2)-0.2).^2 + 0.8 -0.7*x(:,1); %objective function handle
%constraint function handles
Data.constraint{1}=@(x) -exp(x(:,2)-0.2)-x(:,3);
Data.constraint{2}=@(x) x(:,3)+1.1*x(:,1)+1;
Data.constraint{3}=@(x) x(:,2)-1.2*x(:,1);


end%function