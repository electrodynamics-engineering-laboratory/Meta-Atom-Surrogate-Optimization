function rmse=RMSE_calc(Yo,Yp)
%calculates the root mean squared errors between surrogate response and true
%function value
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
%Input:
%Yo - observed function values
%Yp -predicted function values
%Output:
%rmse - root mean squared error
%--------------------------------------------------------------------------


rmse=zeros(size(Yp,2),1);
for ii = 1:size(Yp,2)
   rmse(ii)=sqrt(sum((Yo-Yp(:,ii)).^2));
end
end %function