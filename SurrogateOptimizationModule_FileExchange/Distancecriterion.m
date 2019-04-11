function [ScaledDist,Dist] = Distancecriterion(S,CandPoint)
%calculates distance criterion:
%for every candidate point determine distance to set of already sampled points
%low distance get large scoring value (bad), large distances get small
%scaling value (good)
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
%S - matrix with already sampled points
%CandPoint - matrix with candidate points
%
%Output:
%ScaledDist - distances of candidate points to set of already sampled
%points, scaled to interval [0,1]
%Dist - distances of candidate points to set of already sampled points
%--------------------------------------------------------------------------

%determine distance of each candidate point to set of already sampled
%points
[~,Dist]=knnsearch(S, CandPoint);
MaxDist = max(Dist); %maximum distance of all candidate points to the set of sample points
MinDist = min(Dist); %minimum distance of all candidate points to the set of sample points

if MaxDist ==MinDist
    ScaledDist =ones(length(Dist),1);
else
    ScaledDist = (MaxDist-Dist)/(MaxDist-MinDist); %scaled distances [0,1]
end

end % function