function Data = datainput_hydropower_2plant2_I
%hydropower maximization problem
%two plants, 5 generator types

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
qlow=[750,750,600,600,600]; %lower bounds on capacity of generators
b=[5000;4800];%maximum amounts of water that can be released
Data.objfunction=@(x)poweroutput(x); %objective function handle
%constraint functions
Data.constraint{1}=@(n) n(1)*qlow(1) + n(5)*qlow(5)-b(1);
Data.constraint{2}=@(n) n(2)*qlow(2) + n(3)*qlow(3)+ n(4)*qlow(4)-b(2);
end %function



function y=poweroutput(n) %objective function
qlow=[750,750,600,600,600]; %lower bounds on capacity of generators
qup=[1800,1800,1300,1300,1300]; %upper bounds on capacity of generators
A=[n(1)  0 0 0 n(5);...
   0 n(2) n(3) n(4) 0]; %matrix with generator settings
b=[5000;4800]; %upper bounds on total outflow
x0=qlow+(qup-qlow).*rand(1,5); %initial point for optimizaing water that goes through each generator type
cv=A*qlow(:)-b;%check if with current adjustment there is a feasible solution possible
if any(cv>0) %with the corrent setting of the generators no feasible solution can be generated no matter what the q's are
    y=140+sum(max([zeros(2,1),cv],[],2)); %set objective function value to some positive value
else %otherwise optimize the water that goes through every generator type
    opt=optimset('Display','off');
    [~,y,flag]=fmincon(@(x)objfcn(x,n),x0,A,b,[],[],qlow,qup,@(x)myconstr(x,n),opt);
    if flag<1 %if optimization was not successful, set objective function to some high value
        y=1000;
    end
end
end %poweroutput


function y=objfcn(x,n) %objective function
    %parameter settings for efficiency of single generator types
    c1=[0.5;0.3;0.1;0.035;0.7];
    c2=[-200;10;200;260;-500];
    a=c1.*n';
    b=c2.*n';
    y=a'*x'+b'*ones(5,1);
    y=-y; %maximization problem
end %objfcn

function [c, ceq] =myconstr(x,n) %constraints on total outflow
    c(1,1)=n(1)*x(1)+n(5)*x(5)-5000;
    c(2,1)=n(2)*x(2)+n(3)*x(3) +n(4)*x(4)-4800;    
    ceq=[];
end %myconstr