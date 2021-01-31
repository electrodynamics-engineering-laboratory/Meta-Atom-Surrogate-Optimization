function Data = SOI_OP1(Data,maxeval,P,stdev_int, Surrogate)

%optimization phase 1 for constrained integer optimization problems
%minimize constraint violation function using the candidate point sampling
%approach
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
%Data - structure, contains all problem information
%maxeval - maximum number of allowed function evaluations
%P - perturbation probability for creating candidates
%stdev_int - vector of integers, perturbation ranges
%Surrogate - string, contains the name of the (mixture) surrogate model to
%be used for predictions
%
%Output:
%Data - updated structure with all problem information
%--------------------------------------------------------------------------

[Fbest,idx]=min(Data.pen); % lowest constraint violation value
Data.Fbest_inf = Data.Y(idx); %function value of best infeasible point
Data.xbest_inf=Data.S(idx,:); % best infeasible point 
Data.xbest=Data.xbest_inf;  %best solution so far

% algorithm parameters
NCandidates=1000*Data.dim;       %number of candidate points in each group    
valueweight=0.9; %weight for response surface criterion
distweight=0.1; %weight for distance criterion

%comput constraint violation function values, needed for fitting response
%surface
g_loc=Data.g; 
gpos=sum(max(g_loc,zeros(size(g_loc))),2); %sum up positive constraint violations
gneg=sum(min(g_loc,zeros(size(g_loc))),2);
g_for_surf=gneg;
g_for_surf(gpos>0)=gpos(gpos>0);
Data.Ymed = g_for_surf; %variable Data.Ymed is used to compute surrogate model parameters

%compute response surface parameters using (Data.S, Data.Ymed)
[lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate);

% iterate until max. number of allowed function evaluations reached, or 
%feasible point found
while (size(Data.S,1) < maxeval) && ~Data.Feasible   
    CreateNewCandidates = true;
    while CreateNewCandidates 
        %candidate point generation
        CandPoint = Perturbation_SOI(Data,NCandidates,stdev_int,P);
        %predict function value at candidate points
        CandValue=PredictFunctionValues(Data,Surrogate,CandPoint,lambda,gamma,...
            dmodel,beta,mmodel,w_m);
        %scale predicted objective function values to [0,1]
        ScaledCandValue = PredictedValueCrit(CandValue);
        %compute distance criterion value for each candidate
        [ScaledDist,Dist] = Distancecriterion(Data.S,CandPoint);
        %compute weighted score for every candidate point  
        CandTotalValue = valueweight*ScaledCandValue + distweight*ScaledDist;
        CandTotalValue(Dist <1)=inf; %give bad scores to points that are too close to already sampled points
        if all (CandTotalValue==inf) %regenerate candidate points if necessary
            CreateNewCandidates= true;
        else 
            CreateNewCandidates= false;
        end
    end
    %determine best candidate point's index
    [~,selindex] = min(CandTotalValue);
    % select best candidate point for function evaluation  
    xselected = CandPoint(selindex,:);    
    % clear unnecessary variables
    clear CandPoint CandValue ScaledCandValue Dist ScaledDist CandTotalValue;
    
    % perform function and constraint evaluation at the selected point
    allg=zeros(1,length(Data.constraint));
    gevalst=tic; %record time of constraint evaluations
    for ii = 1:length(Data.constraint)
        allg(:,ii)=feval(Data.constraint{ii},xselected);
    end
    Data.gevaltime=[Data.gevaltime;toc(gevalst)];
    Data.g=[Data.g;allg]; %update matrix with constraint function values
    %compute constraint violation function value by adding up positive
    %constraint function values
    Fselected = sum(max(allg,zeros(1,length(Data.constraint))).^2); 
 
    if  (Fselected < Fbest) %check if improved solution has been found
        Data.xbest_inf = xselected; %update best point found so far
        Data.xbest=Data.xbest_inf;
        Fbest = Fselected; % update best function value found so far
    end
    if Fselected<=0 %if new point is feasible
        Data.Feasible=true; %update feasibility indicator
        Data.xbest=xselected;
    else
        Data.Feasible=false;
    end           
    
    Data.S = [Data.S; xselected]; %update matrix with sample sites
    Data.pen=[Data.pen;Fselected]; %update constraint violation vector
    
    fevalst=tic; %record time for calculationg objective function value
    Data.Y = [Data.Y;feval(Data.objfunction,xselected)]; % update objective function value vector
    Data.fevaltime=[Data.fevaltime;toc(fevalst)];
    
    if Data.Feasible
        Data.fbest=Data.Y(end); %update best function value found so far
        Data.xbest=Data.S(end,:);  %update best point found so far
        break
    else
        %add up positive constraint violations
        g_loc=Data.g;
        gpos=sum(max(g_loc,zeros(size(g_loc))),2);
        gneg=sum(min(g_loc,zeros(size(g_loc))),2);
        g_for_surf=gneg;
        g_for_surf(gpos>0)=gpos(gpos>0);
        Data.Ymed = g_for_surf; %variable Data.Ymed is used to compute surrogate model parameters
        % update response surface model
        [lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate);
    end
    save('Results','Data');
    fprintf('Number of function evaluation: %4.0f; Best infeasible function value: %f\n', size(Data.S,1),Fbest)
end

end %function