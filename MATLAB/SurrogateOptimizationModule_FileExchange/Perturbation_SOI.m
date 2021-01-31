function CandPoint = Perturbation_SOI(Data,NCandidates,stdev_int,P)

%function that generates candidate points for purely integer problems by
%perturbing best point found so far, and by uniformly selecting points from
%the variable domain
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
% Data - structure, contains information about optimization problem
% NCandidates - integer, number of candidate points contained in each group
% stdev_int - integer vectore, perturbation ranges
% P - perturbation probability
%
%Output:
%CandPoint - matrix with candidate points
%--------------------------------------------------------------------------

%set candidate points same as xbest, and then perturb only integer
%variables
CandPoint=repmat(Data.xbest,NCandidates,1);
for ii = 1:NCandidates %for each canidate point
    for jj=1:Data.dim %for each integer variables
        if rand(1)<= P %check if perturbation wanted
            r=randperm(length(stdev_int));  %select random perturbation range
            p=randn(1);
            sign_p=sign(p);
            if sign_p <0
                perturbation= min(-1,round(stdev_int(r(1))*randn(1)));
            else
                perturbation= max(1,round(stdev_int(r(1))*randn(1)));
            end
            CandPoint(ii,jj) = max(Data.xlow(jj),...
                min(CandPoint(ii,jj)+perturbation,...
                Data.xup(jj)));
        end
    end
end

%uniformly selected integer points 
CandPointU= round(repmat(Data.xlow,NCandidates,1) + rand(NCandidates,Data.dim).*repmat(Data.xup-Data.xlow,NCandidates,1));
CandPoint=[CandPoint;CandPointU];

end% function