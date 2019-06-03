function Data = CandidatePointSampling(Data,maxeval,Surrogate,lambda,gamma,...
    dmodel,mmodel,beta, w_m, tolerance)

%optimization by candidate point sampling strategy
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
%Input
%Data - structure, contains all information about the optimization problem
%maxeval - integer, maximum number of allowed function evaluations
%Surrogate - string, surrogate model type to be used
%lambda,gamma  - vectors, parameters of RBF model; lambda=gamma=[] if RBF model not used
%dmodel - structure, parameters for kriging model; dmodel=[] if kriging model not used
%mmodel - structure, parameters for MARS model; mmodel=[] if MARS model not used
%beta - vector, parameters of polynomial regression model; beta=[] if polynomial model not used
%w_m - vector, contains the weights for the models in mixtures; w_m=[] if no mixture model used
%tolerance - scalar, distance when two points are considered equal
%
%Output
%Data - structure, contains updated information about the optimization problem
%--------------------------------------------------------------------------

xrange = Data.xup-Data.xlow; %range between lower and upper variable bounds
minxrange = min(xrange);%smallest variable range
%% perturbation parameters
sigma_stdev(1) = 0.1*minxrange;   %rates for perturbing variables
sigma_stdev(2) = 0.01*minxrange;              
sigma_stdev(3) = 0.001*minxrange;           
NCandPoint=1000*Data.dim;       %number of candidate points in each group    
valueweight=-0.1; %weight for predicted objective function value criterion
%% permutaion probability
if Data.dim>5 %if problem dimension is larger than 5, perturb each variable with prob. 5/dimension
    P=max(0.1,5/Data.dim);
else %if problem dimension is <= 5, perturb all variables with prob. 1
    P=1;
end

iterctr=0;  %initialize iteration counter, needed for determining cycle number   
%% iteration loop
while size(Data.S,1) < maxeval %stopping criterion is max. number allowed function evaluations
    iterctr=iterctr+1;
    % determine which weights to use    
    valueweight=valueweight+.1; %weight for predicted objective function value criterion
    if valueweight>1
        valueweight=0;
    end
    mindistweight=1-valueweight; %weight for distance criterion
    CreateNewCandidates = true;
    while CreateNewCandidates
        %generate candidate points          
        CandPoint=Perturbation(Data, NCandPoint, sigma_stdev,P); 
        %predict objective function values
        CandValue=PredictFunctionValues(Data,Surrogate,CandPoint,lambda,gamma,...
            dmodel,beta,mmodel,w_m);

        % compute scoring criteria to be used for candidate points
        % distances to nearest already sampled point scaled to [0,1]
        % (values close to 0 good)
        [ScaledDist,Dist] = Distancecriterion(Data.S,CandPoint);
        % predicted function values scaled to [0,1] (values close to 0
        % good)
        ScaledCandVal=PredictedValueCrit(CandValue);
        %total weighted score for each candidate point    
        CandTotalValue=valueweight*ScaledCandVal+ mindistweight*ScaledDist;
        % bad scores for points that are too close to already sampled points
        CandTotalValue(Dist < tolerance)=inf;    
        if all (CandTotalValue==inf)
            CreateNewCandidates= true;
        else 
            CreateNewCandidates= false;
        end
    end

    [~,selindex] = min(CandTotalValue); %find index of best candidate point
    %replace MinCandTotalValue by ~ on newer matlab
    xselected = CandPoint(selindex,:); %best candidate point, used for expensive function evaluation        

    % clear unnecessary variables
    clear CandPoint CandValue  CandTotalValue ScaledDist Dist ScaledCandVal;

    % perform function evaluation at the selected point
    fevalst=tic;
    Fselected = feval(Data.objfunction,xselected);
    Data.fevaltime=[Data.fevaltime;toc(fevalst)]; %record time needed for function evaluation

    if (Fselected < Data.fbest) %new point is better than best point found so far
        Data.xbest = xselected; %update best point found so far
        Data.fbest = Fselected; %update best function value found so far
    end

    %update data vectors
    Data.S = [Data.S; xselected]; %sample site matrix
    Data.Y = [Data.Y; Fselected]; %objective function values
    Data.Ymed = Data.Y; %data where large function values set to median, for calculation of surrogate model parameters
    % replace large function values by the median of all available
    % function values
    MedY=median(Data.Y);
    Data.Ymed(Data.Y>MedY)=MedY;

    %recompute model parameters
    [lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate);
    fprintf('Number of function evaluation: %4.0f; Best feasible function value: %f\n', size(Data.S,1),Data.fbest)
    save('Results','Data');
end

end % function