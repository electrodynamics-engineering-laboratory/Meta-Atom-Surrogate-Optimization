function Data= datainput_linearproblem_MI
%replace = by <=0 MINLPLIB
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

Data.xlow=zeros(1,5); %variable lower bounds
Data.xup=[10,10,10,1,1]; %variable upper bounds
%objective function handle
Data.objfunction=@(x)-(- 2*x(:,1) - 3*x(:,2) - 1.5*x(:,3) - 2*x(:,4) + 0.5*x(:,5));
%constraint function handles
Data.constraint{1}=@(x) x(:,1) + x(:,3) - 1.6;
Data.constraint{2}=@(x) 1.333*x(:,2) + x(:,4) - 3;
Data.constraint{3}=@(x) - x(:,3) - x(:,4) + x(:,5); 
Data.integer=(1:3); %indices of integer variables
Data.continuous=(4:5); %indices of continuous variables
Data.dim=5; %problem dimension

end %function