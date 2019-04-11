function Data = datainput_hydropower_1plant3_I
%hydropower maximization problem
%one plant, 5 generator types

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
Data.xup=10*ones(1,5); %variable upper bounds
Data.dim=5; %problem dimension
Data.integer=(1:5); %indices of integer variables
Data.continuous=[]; %indices of continuous variables
qlow=750*ones(1,5); %lower bounds on capacity of generators
b=8000;%maximum amount of water that can be released
Data.objfunction=@(x)poweroutput(x); %objective function handle
Data.constraint{1}=@(n) n(:)'*qlow(:)-b;  %constraint
end %function



function y=poweroutput(n) %objective function
qlow=750*ones(1,5); %lower bounds on capacity of generators
qup=1800*ones(1,5); %upper bounds on capacity of generators
A=n(:)'; %vector with generator settings
b=8000;  %upper bound on total outflow
x0=qlow+(qup-qlow).*rand(1,5); %initial point for optimizaing water that goes through each generator type
cv=A*qlow(:)-b; %check if with current adjustment there is a feasible solution possible
if any(cv>0) %with the corrent setting of the generators no feasible solution can be generated no matter what the q's are
    y=140+(max(0,cv)); %set objective function value to some positive value
else  %otherwise optimize the water that goes through every generator type
    opt=optimset('Display','off');
    [~,y,flag]=fmincon(@(x)objfcn(x,n),x0,A,b,[],[],qlow,qup,@(x)myconstr(x,n),opt);
    if flag<1 %if optimization was not successful, set objective function to some high value
        y=1000;
    end
end
end %poweroutput


function y=objfcn(x,n)%objective function
    %parameter settings for efficiency of single generator types
    c1=[0.6;0.4;0.15;0.25;0.8];
    c2=[-210;50;290;200;-550];
    a=c1.*n';
    b=c2.*n';
    y=a'*x'+b'*ones(5,1);
    y=-y; %maximization problem
end %objfcn

function [c, ceq] =myconstr(x,n) %constraints on total outflow
    c=n(:)*x-8000;
    ceq=[];
end %myconstr