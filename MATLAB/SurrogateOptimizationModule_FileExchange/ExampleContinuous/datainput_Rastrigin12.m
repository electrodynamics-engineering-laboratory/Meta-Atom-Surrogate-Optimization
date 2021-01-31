function Data= datainput_Rastrigin12
%12-dimensional Rastrigin function
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

Data.xlow=-ones(1,12); %lower variable bounds
Data.xup=3*ones(1,12); %upper variable bounds
Data.objfunction=@(x)myfun(x);  %objective function handle
Data.dim = 12;
Data.integer = [];
Data.continuous = (1:12);
end %function

function y=myfun(x)
x=x(:)';
y=sum((x.^2) - cos(2*pi*x),2);
end %myfun