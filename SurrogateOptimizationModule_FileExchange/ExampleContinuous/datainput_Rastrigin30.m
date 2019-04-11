function Data= datainput_Rastrigin30
%30-dimensional Rastrigin function

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

Data.xlow=-ones(1,30); %variable lower bounds
Data.xup=3*ones(1,30); %variable upper bounds
Data.objfunction=@(x)myfun(x); %objective funtion handle
Data.dim = 30; %problem dimension
Data.integer=[];
Data.continuous = (1:30);
end %function

function y=myfun(x) %objective function
x=x(:)';
y=sum((x.^2) - cos(2*pi*x),2);
end %end myfun