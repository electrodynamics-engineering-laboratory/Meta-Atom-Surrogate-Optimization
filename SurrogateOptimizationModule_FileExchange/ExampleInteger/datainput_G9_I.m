function Data=datainput_G9_I
%test problem G9 by Koziel and Michalewicz [1999], integrality constraints
%added
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
Data.objfunction=@(x) (x(:,1)-10).^2+5*(x(:,2)-12).^2 +x(:,3).^4+3*(x(:,4)-11).^2 +...
    10*x(:,5).^6+7*x(:,6).^2+x(:,7).^4-4*x(:,6).*x(:,7)-10*x(:,6)-8*x(:,7);
%define constraints
Data.constraint{1}= @(x) 2*x(:,1).^2+3*x(:,2).^4+x(:,3)+4*x(:,4).^2+5*x(:,5)-127;
Data.constraint{2}= @(x) 7*x(:,1)+3*x(:,2)+10*x(:,3).^2+x(:,4)-x(:,5)-282;
Data.constraint{3}= @(x) 23*x(:,1)+x(:,2).^2+6*x(:,6).^2-8*x(:,7)-196;
Data.constraint{4}= @(x) 4*x(:,1).^2+x(:,2).^2-3*x(:,1).*x(:,2)+2*x(:,3).^2+5*x(:,6)-11*x(:,7);
Data.xlow=-10*ones(1,7); %variable lower bounds
Data.xup=10*ones(1,7); %%variable upper bounds
Data.dim = 7; %problem dimension
Data.integer =(1:7); %indices of integer variables
Data.continuous =[]; %indices of continuous variables 
end %function