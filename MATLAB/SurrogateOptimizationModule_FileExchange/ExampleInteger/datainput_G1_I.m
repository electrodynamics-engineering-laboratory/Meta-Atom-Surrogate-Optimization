function Data=datainput_G1_I
%test problem G1 by Koziel and Michalewicz [1999], integrality constraints
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
Data.objfunction=@(x) 5*sum(x(:,1:4),2) -5*sum(x(:,1:4).^2,2)-sum(x(:,5:13),2);
%define constraints
Data.constraint{1}=@(x)2*x(:,1)+2*x(:,2)+x(:,10)+x(:,11)-10;
Data.constraint{2}=@(x)2*x(:,1)+2*x(:,3)+x(:,10)+x(:,12)-10;
Data.constraint{3}=@(x)2*x(:,2)+2*x(:,3)+x(:,11)+x(:,12)-10;
Data.constraint{4}=@(x)-8*x(:,1)+x(:,10);
Data.constraint{5}=@(x)-8*x(:,2)+x(:,11);
Data.constraint{6}=@(x)-8*x(:,3)+x(:,12);
Data.constraint{7}=@(x)-2*x(:,4)-x(:,5)+x(:,10);
Data.constraint{8}=@(x)-2*x(:,6)-x(:,7)+x(:,11);
Data.constraint{9}=@(x)-2*x(:,8)-x(:,9)+x(:,12);
Data.xlow=zeros(1,13); %variable lower bounds
Data.xup=[ones(1,9),100,100,100,1]; %variable upper bounds
Data.dim = 13; %problem dimension
Data.integer=(1:13); %indices of integer variables
Data.continuous = []; %indices of continuous variables
end %function