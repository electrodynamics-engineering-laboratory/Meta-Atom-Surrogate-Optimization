function [xbest, Fbest]=ODDS(Data, x0,objective, tolerance)

% DDS algorithm for finding the minimum of the response surface...........
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
%input: Data - structure with all problem information
% d - problem dimension
% xbest - starting guess (best point so far)
% objective - function handle to response surface
% tolerance - tolerance when 2 points are considered equal
%output: xbest - best point found during ODDS
%Fbest - response surface value of best point found
%--------------------------------------------------------------------------

xbest=x0;
maxDDSfeval=1000*Data.dim; %maximal number of allowed response surface evaluations
r=0.2; %neighborhood perturbation size parameter

%predicted objective function value at starting point
Fbest=feval(objective, xbest);
Sigma=r*(Data.xup-Data.xlow);  %for perturbation

MaxNoImprovement=20;%maximum number of consecutively failed improvement trials
NoImprovementCounter=0;
fevalctr=1;

while fevalctr<maxDDSfeval && NoImprovementCounter < MaxNoImprovement
    P=1-log(fevalctr)/log(maxDDSfeval);  %perturbation probability
    N=rand(1,Data.dim); %for determining which variables to perturb
    varperturb=find(N<=P); %variables to be perturbed
    if isempty(varperturb) %no variable to be perturbed, enforce one variable to perturb
        c=randperm(Data.dim);  
        varperturb=c(1);
    end
    
    %perturb each variable that is supposed to be perturbed according to
    %DDS algorithm by Tolson and Shoemaker
    for ii = 1:length(varperturb)
        xnew=xbest;
        xnew(varperturb(ii))=xbest(varperturb(ii))+Sigma(varperturb(ii))*randn(1);
        if xnew(varperturb(ii)) <Data.xlow(varperturb(ii))
            xnew(varperturb(ii))=2*Data.xlow(varperturb(ii))-xnew(varperturb(ii));
            if xnew(varperturb(ii)) > Data.xup(varperturb(ii))
                xnew(varperturb(ii))=Data.xlow(varperturb(ii));
            end
        elseif xnew(varperturb(ii))>Data.xup(varperturb(ii))
            xnew(varperturb(ii))=2*Data.xup(varperturb(ii))-xnew(varperturb(ii));
            if xnew(varperturb(ii))<Data.xlow(varperturb(ii))
                xnew(varperturb(ii))=Data.xup(varperturb(ii));
            end
        end
    end
    %check if new point is already contained in set of already evaluated
    %points
    if min(sqrt(sum((repmat(xnew,size(Data.S,1),1)-Data.S).^2,2))) < tolerance
        fnew=inf;
    else %add new point to "already evaluated points" (note that Data structure is NOT an output!!
        fnew=feval(objective,xnew);
        Data.S=[Data.S;xnew];
        Data.Y=[Data.Y;fnew];
    end
    
    if fnew < Fbest %current solution has been improved
        Fbest=fnew; %update best predicted function value
        xbest=xnew; %update best point found
        NoImprovementCounter=0;
    else
        NoImprovementCounter=NoImprovementCounter+1;
    end
    fevalctr=fevalctr+1; %update counter
 
end %while
end%function