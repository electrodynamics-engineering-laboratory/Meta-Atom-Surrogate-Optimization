function Data = datainput_Shekel7

% Shekel-7 function
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

Data.xlow=[0,0,0,0]; %lower variable bounds
Data.xup=[10,10,10,10]; %upper variable bounds
Data.dim=4;
A=[4.0*ones(4,1), ones(4,1), 8*ones(4,1), 6*ones(4,1), [3,2,5,8,6,7;7,9,5,1,2,3.6; ...
    3,2,3,8,6,7;7,9,3,1,2,3.6]];
c=1/10*[1,2,2,4,4,6,3,7,5,5]; 
A=A(:,1:7);
c=c(1:7);
Data.objfunction=@(x)shekel(x,A,c); %handle to objective function
Data.integer=[]; %indices of integer variables
Data.continuous=(1:4); %indices of continuous variables

end %function

function y=shekel(x,A,c)
    x=x(:)';
    y=zeros(size(x,1),1);
    S1=zeros(size(A,2),1);
    for ii =1:size(x,1)
        for jj = 1:size(A,2)
            S1(jj,1)=1./(sum((x(ii,:)'-A(:,jj)).^2)+c(jj));
        end
        y(ii,1)=-sum(S1);
        S1=zeros(size(A,2),1);
    end
end %shekel