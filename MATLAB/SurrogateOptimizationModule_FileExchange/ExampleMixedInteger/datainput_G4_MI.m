function Data=datainput_G4_MI

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
%objective function handle
Data.objfunction=@(x) 5.3578547*x(:,3).^2 +0.8356891*x(:,1).*x(:,5)+37.293239*x(:,1)-40792.141;
%constraint function handles
Data.constraint{1}=@(x) 85.334407+0.0056858*x(:,2).*x(:,5) +0.0006262*x(:,1).*x(:,4)-0.0022053*x(:,3).*x(:,5)-92;
Data.constraint{2}=@(x)-(85.334407+0.0056858*x(:,2).*x(:,5) +0.0006262*x(:,1).*x(:,4)-0.0022053*x(:,3).*x(:,5));
Data.constraint{3}=@(x) 80.51249+0.0071317*x(:,2).*x(:,5)+0.0029955*x(:,1).*x(:,2)+0.0021813*x(:,3).^2-110;
Data.constraint{4}=@(x)-(80.51249+0.0071317*x(:,2).*x(:,5)+0.0029955*x(:,1).*x(:,2)+0.0021813*x(:,3).^2)+90;
Data.constraint{5}=@(x) 9.300961+0.0047026*x(:,3).*x(:,5)+0.0012547*x(:,1).*x(:,3)+0.0019085*x(:,3).*x(:,4)-25;
Data.constraint{6}=@(x)-( 9.300961+0.0047026*x(:,3).*x(:,5)+0.0012547*x(:,1).*x(:,3)+0.0019085*x(:,3).*x(:,4))+20;
Data.xlow=[78,33,27,27,27]; %variable lower bounds
Data.xup=[102,45,45,45,45]; %variable upper bounds
Data.integer=(1:2); %indices of integer variables
Data.continuous=(3:5); %indices of continuous variables
Data.dim = 5; %problem dimension
end %function