function corcoef=cc_calc(Yo,Yp)
%calculates correlation coefficiants of model response and true function
%value
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
%Yp - predicted function values
%Output:
%corcoef - correlation coefficients
%--------------------------------------------------------------------------

corcoef=zeros(size(Yp,2),1);
for ii = 1:size(Yp,2)
    CCmatrix=corrcoef(Yo, Yp(:,ii));
    corcoef(ii)=CCmatrix(1,2);
end
end%function