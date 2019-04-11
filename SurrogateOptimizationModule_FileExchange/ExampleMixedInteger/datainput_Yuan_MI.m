function Data = datainput_Yuan_MI
% test function by Yuan 1988
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
Data.xup=[ones(1,8), 10,10,10]; %variable upper bounds
Data.dim=11;
Data.integer=(1:4); %indices of integer variables
Data.continuous=(5:11); %indices of continuous variables
%objective function handle
Data.objfunction=@(x) (x(:,5)-1).^2 + (x(:,6)-2).^2 + (x(:,7)-1).^2 -...
    log(x(:,8)+1) + (x(:,9)-1).^2 + (x(:,10)-2).^2 + (x(:,11)-3).^2;
%constraintfunction handles
Data.constraint{1}=@(x) x(:,1) +x(:,2)+x(:,3)+x(:,9)+x(:,10)+x(:,11)-5;
Data.constraint{2}=@(x) x(:,7).^2+x(:,9).^2+x(:,10).^2+x(:,11).^2-5.5;
Data.constraint{3}=@(x) x(:,1)+x(:,9)-1.2;
Data.constraint{4}=@(x) x(:,2)+x(:,10)-1.8;
Data.constraint{5}=@(x) x(:,3)+x(:,11)-2.5;
Data.constraint{6}=@(x) x(:,4)+x(:,9)-1.2;
Data.constraint{7}=@(x) x(:,6).^2+x(:,10).^2-1.64;
Data.constraint{8}=@(x) x(:,7).^2+x(:,11).^2-4.25;
Data.constraint{9}=@(x) x(:,6).^2+x(:,11).^2-4.64;
Data.constraint{10}=@(x) x(:,5)-x(:,1);
Data.constraint{11}=@(x) x(:,6)-x(:,2);
Data.constraint{12}=@(x) x(:,7)-x(:,3);
Data.constraint{13}=@(x) x(:,8)-x(:,4);

end %function 