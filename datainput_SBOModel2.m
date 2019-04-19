function Data = datainput_SOBModelProject2
%Voltage divider problem with polar coordinates

Data.xlow=[0,0,0,0]; %variable lower bounds in cartesian to be converted to polar
Data.xup=[10,10,10,10]; %variable upper bounds in cartesian to be converted to polar
%objective function
Data.objfunction= @(x)(sqrt(x(:,3)^2+x(:,4)^2)*exp(1i.*atan2(x(:,4),x(:,3)))./(sqrt(x(:,1)^2+x(:,2)^2)*exp(1i.*atan2(x(:,2),x(:,1))) + sqrt(x(:,3)^2+x(:,4)^2)*exp(1i.*atan2(x(:,4),x(:,3)))));
Data.integer=[]; %indices of integer variables
Data.continuous=(1:4); %indices of continuous variables
Data.dim = 4; %problem dimension
end %function
