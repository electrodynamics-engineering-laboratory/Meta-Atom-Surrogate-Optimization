function xbest=Multi_DDS(Data, valueweight,mindistweight, tolerance, myfun)

%%DDS algorithm that generates from one point N perturbations, computes the
%%weighted criterion for each of these points, and selects  the best from 
%these N points as starting point for next perturbations
%
%initial point is the best point found so far
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
%Data - structure that contains all information about the problem
%valueweight,mindistweight - weights for scoring criteria
%tolerance - distance at which 2 points are considered equal
%myfun - handle to the (mixture) surrogate model
%
%Output:
%xbest - the best point found during the optimization
%--------------------------------------------------------------------------

maxfeval=Data.dim*1000; %maximum number of allowed response surface evaluations
NumberOffspring= 20;%dim*10; %number of points created by DDS perturbation in each iteration
MaxIter=maxfeval/NumberOffspring; %max number of iterations allowed

r=0.2; %neighborhood perturbation size parameter
Sigma=r*(Data.xup-Data.xlow);
iterctr=0;
old_offspring=[];
old_prediction=[];
no_improved_child =0;
xbest=Data.xbest;

while iterctr<MaxIter
    iterctr=iterctr+1;
    P=1-log(iterctr)/log(MaxIter);  %perturbation probability for each variable
    offspring=zeros(NumberOffspring,Data.dim); %initialize array for offspring
    Fpred=zeros(NumberOffspring,1); %initialize aray for predicted function values
    for jj= 1:NumberOffspring    
        N=rand(1,Data.dim); %create random vector, for deciding which variables to perturb
        varperturb=find(N<=P); %variables to be perturbed
        if isempty(varperturb) %if no variable is supposed to be perturbed, force that one randomly picked
                                %variable will be perturbed
            c=randperm(Data.dim);
            varperturb=c(1);
        end
        xnew=Data.xbest; %initialize new point as best point
        %do perturbations according to DDS algorithm
        for ii = 1:length(varperturb)            
            xnew(varperturb(ii))=xnew(varperturb(ii))+Sigma(varperturb(ii))*randn(1);
            if xnew(varperturb(ii)) <Data.xlow(varperturb(ii)) %reflect perturbation
                xnew(varperturb(ii))=2*Data.xlow(varperturb(ii))-xnew(varperturb(ii));
                if xnew(varperturb(ii)) > Data.xup(varperturb(ii))
                    xnew(varperturb(ii))=Data.xlow(varperturb(ii));
                end
            elseif xnew(varperturb(ii))>Data.xup(varperturb(ii)) %reflect perturbation
                xnew(varperturb(ii))=2*Data.xup(varperturb(ii))-xnew(varperturb(ii));
                if xnew(varperturb(ii))<Data.xlow(varperturb(ii))
                    xnew(varperturb(ii))=Data.xup(varperturb(ii));
                end
            end
        end
        offspring(jj,:)=xnew;
        Fpred(jj,1)=myfun(offspring(jj,:));
    end
    old_offspring=[old_offspring;offspring];
    old_prediction=[old_prediction;Fpred];
    %for every offspring predict objective function value with response
    %surface, and compute distance to set of already sampled points
    [ScaledCandMinDist,CandMinDist] = Distancecriterion(Data.S,old_offspring);
    
    % predicted function values scaled to [0,1], values close to 0
    % are good
    ScaledCandVal=PredictedValueCrit(old_prediction);
    %total weighted score for each candidate point    
    CandTotalValue=valueweight*ScaledCandVal+ mindistweight*ScaledCandMinDist;
    % bad scores for points that are too close to already sampled points
    CandTotalValue(CandMinDist < tolerance)=inf;
    
    %CandTotalValueIter=CandTotalValue(end-NumberOffspring+1:end); %these are the scores for the current population
    [~,id]=min(CandTotalValue);
    xbest_old=xbest;
    xbest=old_offspring(id, :);
    
    if (iterctr-1)*NumberOffspring >=id  || sqrt(sum((xbest-xbest_old).^2,2))< tolerance%&& id <=(iterctr+1)*NumberOffspring
        no_improved_child =no_improved_child +1; %current best point could not be improved
    else
        no_improved_child =0; %current best point could be improved
    end
    if no_improved_child >5 %stop if 5 consecutive trials did not give improved point
        break
    end     
end

end%function