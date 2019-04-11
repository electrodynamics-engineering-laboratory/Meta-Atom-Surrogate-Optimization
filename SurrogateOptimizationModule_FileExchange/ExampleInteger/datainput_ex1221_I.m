function Data= datainput_ex1221_I
%problem ex1221 given on http://www.gamsworld.org/minlp/minlplib/ex1221.htm 
%integrality constraints added, = replaced by <=

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
Data.objfunction=@(x)-(-2*x(:,1) - 3*x(:,2) - 1.5*x(:,3) - 2*x(:,4) + 0.5*x(:,5));
%define constraints
Data.constraint{1}=@(x) x(:,1).^2 + x(:,3) - 1.25;
Data.constraint{2}=@(x) x(:,2).^1.5 + 1.5*x(:,4) -3;
Data.constraint{3}=@(x) x(:,1) + x(:,3) - 1.6;
Data.constraint{4}=@(x) 1.333*x(:,2) + x(:,4) - 3;
Data.constraint{5}=@(x) - x(:,3) - x(:,4) + x(:,5); 
Data.xlow=zeros(1,5); %variable lower bounds
Data.xup=[10,10,10,1,1]; %variable upper bounds
Data.dim=5; %problem dimension
Data.integer = (1:5); %indices of integer variables
Data.continuous  = []; %indices of continuous variables
end %function