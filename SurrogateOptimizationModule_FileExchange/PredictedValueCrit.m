function ScaledCandValue = PredictedValueCrit(CandValue)
%scales predicted objective function values at candidate points to [0,1]
%small values are good, small value gets low number close to 0 in scaling
%
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
%CandValue - vector with predicted objective function values
%
%Output:
%ScaledCandValue - vector with predicted objective function values scaled
%to [0,1]
%--------------------------------------------------------------------------

MinCandValue = min(CandValue); %find smallest predicted value
MaxCandValue = max(CandValue); %find largest predicted value
if MinCandValue ==MaxCandValue
    ScaledCandValue =ones(length(CandValue),1);
else
    ScaledCandValue = (CandValue-MinCandValue)/(MaxCandValue-MinCandValue);
end
end %function