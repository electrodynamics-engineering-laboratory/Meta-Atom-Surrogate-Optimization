function Data = datainput_BermanAshrafi_MI

%problem by Berman and Ashrafi, 1993
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

Data.xlow=zeros(1,11); %variable lower bounds
Data.xup=[ones(1,8),0.997,0.9985,0.9988]; %variable upper bounds
Data.dim=11;
Data.integer=(1:4); %indices of integer variables
Data.continuous=(5:11); %indices of continuous variables

Data.objfunction=@(x) -x(:,9).*x(:,10).*x(:,11); %objective function handle
%constraint function handles
Data.constraint{1}=@(x) -x(:,1)-x(:,2)-x(:,3)+1;
Data.constraint{2}=@(x) -x(:,4)-x(:,5)-x(:,6)+1;
Data.constraint{3}=@(x) -x(:,7)-x(:,8)+1;
Data.constraint{4}=@(x) 3*x(:,1) +x(:,2)+2*x(:,3)+3*x(:,4) +2*x(:,5)+...
    x(:,6)+3*x(:,7) +2*x(:,8)-10;
Data.constraint{5}=@(x) log(0.1)*x(:,1)+log(0.2)*x(:,2) +log(0.15)*x(:,3)-...
    log(1-x(:,9));
Data.constraint{6}=@(x) log(0.05)*x(:,4)+log(0.2)*x(:,5) +log(0.15)*x(:,6)-...
    log(1-x(:,10));
Data.constraint{7}=@(x) log(0.02)*x(:,7)+log(0.06)*x(:,8) -log(1-x(:,11));

end %function