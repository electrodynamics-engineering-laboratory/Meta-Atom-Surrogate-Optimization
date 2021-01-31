function Data = SOI_OP2(Data,maxeval,P,stdev_int,Surrogate)

%optimization phase 2
%try to improve the best feasible point found so far
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
% Data - structure, contains all problem information
% maxeval - maximum number of allowed function evaluations
% P - perturbation probability for creating canididate points
% stdev_int - vector of integers, perturbation ranges
% Surrogate - string, contains name of the (mixture) surrogate model to be
% used
%
%Output:
%Data - updated structure with problem information
%--------------------------------------------------------------------------

%parameters
NCandidates=1000*Data.dim;       %number of candidate points in each group    
valueweight=-0.1; %weight for predicted objective function value criterion

if isfield(Data,'constraint') %constrained problems
    feas=find(Data.pen <=0); %check which points are feasible
    infeas=find(Data.pen >0);%check which points are infeasible
    Sfeas=Data.S(feas,:);    %set of feasible points
    Yfeas=Data.Y(feas);      %function values of feasible points
    [Data.fbest,id]=min(Yfeas);   %find best feasible function value
    Data.xbest=Sfeas(id,:);       %find best feasible point
    Data.Ypenalized=zeros(length(Data.Y),1); %adjust function values of infeasible points
    Data.Ypenalized(feas)=Data.Y(feas);
    Data.Ypenalized(infeas)= max(Yfeas)+100*Data.pen(infeas);
else %unconstrained problems
    [Data.fbest,idx]=min(Data.Y); %find best feasible function value
    Data.xbest=Data.S(idx,:);     %find best feasible point
    Data.Ypenalized=Data.Y;  
end

Data.Ymed=Data.Ypenalized; %vector where large function values are replaced by median
MedY=median(Data.Ypenalized); %median of function values 
Data.Ymed(Data.Ypenalized>MedY)=MedY; %replace large function values by median, used in calculation of surrogate model parameters

%compute parameters of response surface
[lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate);
   
% iterate until max number of function evaluations reached
while size(Data.S,1) < maxeval %stopping criterion is max. number allowed function evaluations
    % determine which weights to use    
    valueweight=valueweight+.1; %weight for predicted objective function value criterion
    if valueweight>1
        valueweight=0;
    end
    mindistweight=1-valueweight; %weight for distance criterion
    
    CreateNewCandidates = true;
    while CreateNewCandidates
        %generate candidate points          
        CandPoint=Perturbation_SOI(Data,NCandidates,stdev_int,P);
        %predict objective function values
        CandValue=PredictFunctionValues(Data,Surrogate,CandPoint,lambda,gamma,...
            dmodel,beta,mmodel,w_m);
        % values for distance criterion
        [ScaledDist,Dist] = Distancecriterion(Data.S,CandPoint);
        % predicted function values scaled to [0,1] 
        ScaledCandVal=PredictedValueCrit(CandValue);
        %total weighted score for each candidate point    
        CandTotalValue=valueweight*ScaledCandVal+ mindistweight*ScaledDist;
        CandTotalValue(Dist < 1)=inf;% bad scores for points that are too close to already sampled points    
        if all (CandTotalValue==inf) %create new candidates if these ones have already been sampled
            CreateNewCandidates= true;
        else 
            CreateNewCandidates= false;
        end
    end
    
    [~,selindex] = min(CandTotalValue); %find index of best candidate point
    xselected = CandPoint(selindex,:); %best candidate point, used for expensive function evaluation        

    % clear unnecessary variables
    clear CandPoint CandValue  CandTotalValue ScaledDist Dist ScaledCandVal;

    % perform function evaluation at the selected point
    fevalst=tic;
    Fselected = feval(Data.objfunction,xselected);
    Data.fevaltime=[Data.fevaltime;toc(fevalst)]; %record time needed for function evaluation
    
    % do constraint evaluations if necessary
    if isfield(Data,'constraint')
        Gselected=zeros(1,length(Data.constraint));   
        gevalst=tic; %record time required for constraint evaluation
        for ii =1:length(Data.constraint)
            Gselected(1,ii)=feval(Data.constraint{ii},xselected);
        end
        Data.gevaltime=[Data.gevaltime;toc(gevalst)];
        Data.g=[Data.g;Gselected];
        Data.pen=[Data.pen; sum(max(zeros(1,length(Data.constraint)),Gselected).^2)];
    end
    
    if (Fselected < Data.fbest)   %check if improvement found
        if isfield(Data,'constraint') && Data.pen(end)<=0 %for constraint problems check feasibility
            Data.xbest = xselected; %update best point found so far
            Data.fbest = Fselected;
        elseif ~isfield(Data,'constraint') %unconstrained problems
            Data.xbest = xselected;%update best point found so far
            Data.fbest = Fselected;
        end
    end
    Data.S = [Data.S; xselected]; %update sample site matrix
    Data.Y = [Data.Y; Fselected]; %update function value vector
    if isfield(Data,'constraint') %for constrained problems
        infeas=find(Data.pen>0);  %find infeasible points
        feas=find(Data.pen <=0);  %find feasible points
    end
    
    if isfield(Data,'constraint') %equip function values with penalty for infeasible points
        Data.Ypenalized=zeros(length(Data.Y),1);
        Data.Ypenalized(feas)=Data.Y(feas);
        Data.Ypenalized(infeas)= max(Yfeas)+100*Data.pen(infeas);
    else
        Data.Ypenalized=Data.Y;
    end
    
    Data.Ymed=Data.Ypenalized; %vector where large function values are replaced by median
    MedY=median(Data.Ypenalized); %median of function values 
    Data.Ymed(Data.Ypenalized>MedY)=MedY; %replace large function values by median, used in calculation of surrogate model parameters

    %update model parameters
    [lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate);
    fprintf('Number of function evaluation: %4.0f; Best feasible function value: %f\n', size(Data.S,1),Data.fbest)
    save('Results','Data');
end

end % function