function Data= datainput_throughput_I
%23-dimensional throughput maximization problem
%Pichitlamken et al [2006]
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

Data.xlow=ones(1,23); %variable lower bounds
Data.xup=20*ones(1,23); %variable upper bounds
%objective function
Data.objfunction=@(x) myobjfunction(x);
%constraints
Data.constraint{1}=@(x) sum(x(:,(1:11)),2)-80;
Data.constraint{2}=@(x) sum(x(:,(12:23)),2)-80;
Data.integer =(1:23); %indices of integer variables
Data.continuous=[]; %indices of continuous variables
Data.dim = 23; %problem dimension
end %function

function y=myobjfunction(x) %objective function
b=x(:,1:11); %buffer sizes
r=x(12:23);  %service rates
%parameters
runlength=10;
seed=2;
f=Throughput(r, runlength, seed, b); %objective function simulation
y=-f; %maximization problem
end %myobjfunction