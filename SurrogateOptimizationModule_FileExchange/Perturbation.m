function CandPoint = Perturbation(Data, NCandPoint,  sigma_stdev,P)

% generate candidate points for next sample site by perturbing vaiables
% obtain two groups of candidate points, each with  1000*dim points
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
% Data - structure with all problem data
% NCandPoint=1000*dim - number of candidate points generated for each group
% sigma_stdev - perturbation rate for continuous variables
% P - probability with which each variable is perturbed
% dim - problem dimension
%Output:
% CandPoint - candidate points
%--------------------------------------------------------------------------

%Group 1: perturbations of best point found so far
CandPoint_G1=repmat(Data.xbest, NCandPoint,1);
for ii = 1:NCandPoint %generate NCandPoint candidate points
    for jj = 1:Data.dim %go through each dimension
        if rand(1) < P %check if perturbation should be done
            p=randperm(length(sigma_stdev)); %randomly pick the rate for pertubation
            CandPoint_G1(ii,jj) = max(Data.xlow(jj), ...
                min(CandPoint_G1(ii,jj)+sigma_stdev(p(1))*randn(1),...
                Data.xup(jj)));
        end
    end
end
    
%Group 2: randomly pick variabe values from whole variable domain
CandPoint_G2=repmat(Data.xlow,NCandPoint,1) + ...
    rand(NCandPoint,Data.dim).* repmat(Data.xup-Data.xlow,NCandPoint,1);
CandPoint=[CandPoint_G1;CandPoint_G2];

clear CandPoint_G1 CandPoint_G2;

end %function