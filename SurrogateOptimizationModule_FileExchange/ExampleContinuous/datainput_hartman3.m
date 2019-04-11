function Data=datainput_hartman3
 
%3-dimensional hartmann function
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
Data.xlow =zeros(1,3); % variable lower bounds
Data.xup=ones(1,3);     % variable upper bounds
Data.objfunction=@(x)myfun(x); %handle to objective function
Data.dim = 3; %problem dimesnion
Data.integer = []; %indices of integer variables
Data.continuous = (1:3); %indices of continuous variables

end %function

function y=myfun(x) %objective function
c=[1,1.2,3,3.2]'; %c,A,b are data vectors
A=[3, 10, 30; 0.1,10,35; 3, 10, 30;0.1,10,35];
P=[0.3689    0.1170    0.2673
    0.4699    0.4387    0.7470
    0.1091    0.8732    0.5547
    0.0382    0.5743    0.8828];
x=x(:)'; % make sure vector is row vector
y=-sum(c.*exp(-sum(A.*(repmat(x,4,1)-P).^2,2))); %compute objective function value
end %myfun