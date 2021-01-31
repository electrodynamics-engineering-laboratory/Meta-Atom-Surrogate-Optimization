function Data = datainput_Powell24
%24-dimesnional Powell function

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

Data.xlow=-4*ones(1,24); %lower variable bounds
Data.xup=5*ones(1,24); %upper variable bounds
Data.dim = 24; %problem dimension
Data.objfunction= @(x) myfun(x); %handle to objective function
Data.integer = []; %indices of integer variables
Data.continuous =(1:24); %indices of continuous variables
end %function

function y=myfun(x)
n = 24;
m = n;
for i = 1:m/4
    fvec(4*i-3) = x(4*i-3)+10*(x(4*i-2));
    fvec(4*i-2) = sqrt(5)*(x(4*i-1)-x(4*i));
    fvec(4*i-1) = (x(4*i-2)-2*(x(4*i-1)))^2;
    fvec(4*i)   = sqrt(10)*(x(4*i-3)-x(4*i))^2;
end;
fvec = fvec';
y = norm(fvec)^2;
end %myfun