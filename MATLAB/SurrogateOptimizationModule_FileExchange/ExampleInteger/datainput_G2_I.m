function Data=datainput_G2_I
% test problem G2 by Koziel and Michalewicz [1999], integrality constraints
% added
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
Data.dim=25; %problem dimension
%define objective function
Data.objfunction = @(x) -abs((sum(cos(x).^4,2)-2*prod(cos(x).^2,2))./(sqrt(sum(repmat((1:Data.dim),size(x,1),1).*x.^2,2))) );
%define constraints
Data.constraint{1}= @(x) -prod(x,2)+0.75;
Data.constraint{2}= @(x) sum(x,2)-7.5*Data.dim;
Data.xlow=zeros(1,Data.dim); %variable lower bounds
Data.xup=10*ones(1,Data.dim); %variable upper bounds
Data.integer = (1:25); %indices of integer variables
Data.continuous =[]; %indices of continuous variables
end %function