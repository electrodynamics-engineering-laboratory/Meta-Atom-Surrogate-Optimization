function Data= datainput_throughput_small_I
%small scale throughput maximization problem
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

Data.xlow=ones(1,5); %variable lower bounds
Data.xup=20*ones(1,5); %variable upper bounds
%objective function
Data.objfunction=@(x) myobjfunction(x);
%constraints
Data.constraint{1}=@(x) sum(x(:,(1:2)),2)-20;
Data.constraint{2}=@(x) sum(x(:,(3:5)),2)-20;
Data.integer =(1:5); %indices of integer variables
Data.continuous=[]; %indices of continuous variables
Data.dim = 5; %problem dimension
end %function

function y=myobjfunction(x) %objective function
b=x(:,1:2); %buffer sizes
r=x(3:5);   %service rates
%parameters
runlength=10; 
seed=2;
f=Throughput(r, runlength, seed, b); %throughput simulation
y=-f; %maximization problem
end %myobjfunction