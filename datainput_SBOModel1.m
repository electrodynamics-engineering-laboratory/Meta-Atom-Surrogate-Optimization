function Data = datainput_SOBModelProject1
%Voltage divider problem input with cartesian inputs

Data.xlow=[0,-10,0,-10]; %variable lower bounds
Data.xup=[10,10,10,10]; %variable upper bounds
%objective function
Data.objfunction= @(x)((x(:,3)+i*x(:,4))./(x(:,1)+i*x(:,2) + x(:,3)+i*x(:,4)));
Data.integer=[]; %indices of integer variables
Data.continuous=(1:4); %indices of continuous variables
Data.dim = 4; %problem dimension
end %function
