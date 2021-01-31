function mae= MAE_cal(Yo,Yp)
%calculates the maximal absolute error between surrogate response and true
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
%
%Input:
%Yo - original objective function value
%Yp - predicted objective function value
%Output:
%mae - maximal absolute error
%--------------------------------------------------------------------------

mae=zeros(size(Yp,2),1);
for ii = 1:size(Yp,2)
    mae(ii)=max(abs(Yo-Yp(:,ii)));
end
end %function