function mad=MAD_cal(Yo,Yp)
%calculates median absolute difference between surrogate response and true
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
%Yo - observed objective function value
%Yp - predicted objective function value
%Output:
%mad - median absolute difference
%--------------------------------------------------------------------------

mad=zeros(size(Yp,2),1);
for ii = 1:size(Yp,2)
    mad(ii)=median(abs(Yo-Yp(:,ii)));
end
end%function