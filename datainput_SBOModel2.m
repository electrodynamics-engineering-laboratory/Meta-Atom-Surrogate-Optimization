function Data = datainput_SOBModelProject
%Voltage divider problem with polar coordinates

Data.xlow=[-10,-10,-10,-10]; %variable lower bounds
Data.xup=[10,10,10,10]; %variable upper bounds
%objective function
Data.objfunction= @(x)(sqrt(x(:,3)^2+x(:,4)^2)*exp(atan2(x(:,3),x(:,4)))./(sqrt(x(:,1)^2+x(:,2)^2)*exp(atan2(x(:,2),x(:,1))) + sqrt(x(:,3)^2+x(:,4)^2)*exp(atan2(x(:,4),x(:,3))));
Data.integer=[]; %indices of integer variables
Data.continuous=(1:4); %indices of continuous variables
Data.dim = 4; %problem dimension
end %function