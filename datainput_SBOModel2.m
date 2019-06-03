function Data = datainput_SOBModelProject2
%Voltage divider problem with polar coordinates

Data.xlow=[0,0,0,0]; %variable lower bounds in cartesian to be converted to polar
Data.xup=[10,2*pi,10,2*pi]; %variable upper bounds in cartesian to be converted to polar
%objective function
Data.objfunction= @(x)(x(:,3)*exp(i*x(:,4))/(x(:,1)*exp(i*x(:,2) + x(:,3)*exp(i*x(:,4)))));
Data.integer=[]; %indices of integer variables
Data.continuous=(1:4); %indices of continuous variables
Data.dim = 4; %problem dimension
end %function
