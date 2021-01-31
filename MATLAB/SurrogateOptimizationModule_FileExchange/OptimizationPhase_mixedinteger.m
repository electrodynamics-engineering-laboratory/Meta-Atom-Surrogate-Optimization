function Data= OptimizationPhase_mixedinteger(Data,Surrogate,maxeval)

%OptimizationPhase_MI iteratively generates new sample points for
%evaluating the costly objective function. In every iteration 4 new sample
%points are evaluated and the response surface updated. The sample points
%are generated as follows (candidate point sampling):
% a) pertub integer variables of best point found so far
% b) pertub continuous variables of best point found so far
% c) pertub integer and continuous variables of best point found so far
% d) uniformly selected points from variable domain
% For each group 500*dimension candidates are generated and a weighted
% scoring criterion is used to select the best point from each group
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
%Data - structure with all problem information 
%Surrogate - string, contains name (mixture) surrogate model to be used
%maxeval - integer, maximum number of allowed function evaluations
%
%Output:
%Data - updated Data structure
%--------------------------------------------------------------------------

%determine perturbation probability for each variable when generating
%candidate points
if Data.dim>5 %if problem dimension is larger than 5, perturb each variable with prob. 5/dimension
    P=max(0.1,5/Data.dim);
else %if problem dimension is <= 5, perturb all variables with prob. 1
    P=1;
end

xrange = Data.xup-Data.xlow; %range between lower and upper variable bounds
minxrange = min(xrange);%smallest variable range
tolerance = 1e-6; %distance tolerance for 2 points to be considered equal

stdev_cont(1) = 0.1*minxrange;   %rates for perturbing continuous variables
stdev_cont(2) = 0.01*minxrange;  %rates for perturbing continuous variables          
stdev_cont(3) = 0.001*minxrange; %rates for perturbing continuous variables
stdev_int(1)=max(1, round(stdev_cont(1))); %rates for perturbing integer variables
stdev_int(2)=max(1, round(stdev_cont(2))); %rates for perturbing integer variables    
stdev_int(3)=max(1, round(stdev_cont(3))); %rates for perturbing integer variables    

NCandidates=500*Data.dim;       %number of candidate points in each group    


%compute constraint function values at points in starting design
if isfield(Data,'constraint') %check if problem is constrained
    g_val=zeros(size(Data.S,1),length(Data.constraint));
    Data.gevaltime=zeros(size(Data.Y));
    for jj=1:size(Data.S,1) %check for every point if constraints satisfied
        gevalst=tic; %record the time of constraint evaluations
        for ii = 1:length(Data.constraint)
            g_val(jj,ii)=feval(Data.constraint{ii},Data.S(jj,:));
        end
        Data.gevaltime(ii,1)=toc(gevalst);
    end
    Data.g=g_val; %matrix containing constraint function values
    % sum up squares of constraint function values that are larger than 0 
    Data.pen=sum(max(zeros(size(g_val)),g_val).^2,2); 
    feas=find(Data.pen <=0); %indices of feasible sample points
    infeas=find(Data.pen >0); %indices of infeasible sample points
    Sfeas=Data.S(feas,:); %feasible sample sites
    Yfeas=Data.Y(feas);   %feasible objective function values
    [Data.fbest,id]=min(Yfeas); %best feasible objective function value
    Data.xbest=Sfeas(id,:);  %best feasible sample site
    Data.Ypenalized=zeros(length(Data.Y),1); %initialize vector with penalized objective function values
    Data.Ypenalized(feas)=Data.Y(feas); 
    Data.Ypenalized(infeas)= max(Yfeas)+100*Data.pen(infeas);
else %there are no constraints
    [Data.fbest,idx]=min(Data.Y);  %best function value so far
    Data.xbest=Data.S(idx,:); %best point found so far
    Data.Ypenalized=Data.Y;  
end
    

Data.Ymed=Data.Ypenalized; %vector where large function values are replaced by median
MedY=median(Data.Ypenalized); %median of function values 
Data.Ymed(Data.Ypenalized>MedY)=MedY; %replace large function values by median, used in calculation of surrogate model parameters

%compute parameters of (mixture) surrogate models
[lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate);

valueweight=-.1;%initialize weight for response surface criterion

% iteration loop
while size(Data.S,1) < maxeval %do until max. number of allowed function evaluations reached       
    %update weights for scoring criteria
    valueweight=valueweight+0.1;
    if valueweight >1
        valueweight=0;
    end
    distweight=1-valueweight; %weight for the distance criterion
    
    %create candidate points until at least one point generated that is not
    %yet sampled
    CreateNewCandidates = true;
    while CreateNewCandidates
        %create candidate points by perturbing discrete, continuous, discrete
        %and continuous variables, and by randomly sampling points over domain
        CandPoint = Perturbation_SOMI(Data,NCandidates,stdev_int,stdev_cont,P);
        %predict objective function values of candidate points
        CandValue=PredictFunctionValues(Data,Surrogate,CandPoint,lambda,gamma,...
                dmodel,beta,mmodel,w_m);
        %compute distance and scaled distances to set of already sampled
        %points
        [ScaledDist,Dist] = Distancecriterion(Data.S,CandPoint);
        %scale predicted objective function values to [0,1]
        ScaledCandValue = PredictedValueCrit(CandValue);   
        % compute total weighted score of criterion 1 and 2, low numbers are
        % good
        CandTotalValue=valueweight*ScaledCandValue+ distweight*ScaledDist;
        % bad scores for points that are too close to already sampled points
        CandTotalValue(Dist < tolerance)=inf;
        if all (CandTotalValue==inf) %regenerate candidate points if necessary
            CreateNewCandidates= true;
        else 
            CreateNewCandidates= false;
        end
    end
        
        
    CandPoint_1=CandPoint(1:500*Data.dim,:); %all points in group 1
    CandTotalValue_1=CandTotalValue(1:500*Data.dim); %all total scores of points in group 1
    [~,selindex(1)] = min(CandTotalValue_1); %find best point in group 1
    if CandTotalValue_1(selindex(1)) == inf
        xselected(1,:)=inf*ones(1,Data.dim);
    else
        xselected(1,:)= CandPoint_1(selindex(1),:); %best point in group 1 
    end
    CandPoint_2=CandPoint(500*Data.dim+1:1000*Data.dim,:); %all points in group 2
    CandTotalValue_2=CandTotalValue(500*Data.dim+1:1000*Data.dim); %all total scores of points in group 2
    [~,selindex(2)] = min(CandTotalValue_2); %find best point in group 2
    if CandTotalValue_2(selindex(2)) == inf
        xselected(2,:)=inf*ones(1,Data.dim);
    else
        xselected(2,:)= CandPoint_2(selindex(2),:); %best point in group 2 
    end
    CandPoint_3=CandPoint(1000*Data.dim+1:1500*Data.dim,:);%all points in group 3
    CandTotalValue_3=CandTotalValue(1000*Data.dim+1:1500*Data.dim);  %all total scores of points in group 3
    [~,selindex(3)] = min(CandTotalValue_3); %find best point in group 3
    if CandTotalValue_3(selindex(3)) == inf
        xselected(3,:)=inf*ones(1,Data.dim);
    else
        xselected(3,:)= CandPoint_3(selindex(3),:); %best point in group 3 
    end
    CandPoint_4=CandPoint(1500*Data.dim+1:2000*Data.dim,:); %all points in group 4
    CandTotalValue_4=CandTotalValue(1500*Data.dim+1:2000*Data.dim); %all total scores of points in group 4
    [~,selindex(4)] = min(CandTotalValue_4); %find best point in group 4
    if CandTotalValue_4(selindex(4)) == inf
        xselected(4,:)=inf*ones(1,Data.dim);
    else
        xselected(4,:)= CandPoint_4(selindex(4),:); %best point in group 4 
    end
    %check if any of the new points is too close to an already sampled
    %points, and delete it if necessary (may hapen whe pertrubing integers
    %only, and the range of integer variables is very small)
    ii =1; %initialize counter
    no_xnew=size(xselected,1); %number of new sample sites
    while ii<=no_xnew 
        if all(xselected(ii,:)==inf) %in case the selected point is too close to already sampled point
            xselected(ii,:)=[]; %delete this point
            no_xnew=no_xnew-1;
        else 
            ii = ii+1;
        end
    end
           
    % clear unnecessary variables
    clear CandPoint CandValue  CandTotalValue ScaledDist Dist...
        ScaledCandValue CandPoint_1 CandPoint_2 CandPoint_3 CandPoint_4...
        CandTotalValue_1 CandTotalValue_2 CandTotalValue_3 CandTotalValue_4
    
    % perform expensive function evaluation at the selected point
    Fselected=zeros(no_xnew,1);
    %sequentially evaluating new points
    for ii=1:no_xnew
        fevalst=tic;
        Fselected(ii,1) = feval(Data.objfunction,xselected(ii,:));
        Data.fevaltime=[Data.fevaltime;toc(fevalst)];
    end
    
    %parallel matlab to evaluate new points %define handle to objf
    %time=zeros(no_xnew,1);
%     parfor ii=1:no_xnew
%         fevalst=tic; %record computation time to obtain funvtion evaluations
%         Fselected(ii,1) = feval(objf,xselected(ii,:)); %expensive function evaluation
%         time(ii)=toc(fevalst); 
%     end
%     Data.fevaltime=[Data.fevaltime;time];
        
    
    if isfield(Data,'constraint')%check constraints at new points
        for jj=1:no_xnew 
            Gselected=zeros(1,length(Data.constraint));   
            gevalst=tic; %record time required for constraint evaluation
            for ii = 1:length(Data.constraint) %compute constraint function values 
                Gselected(1,ii)=feval(Data.constraint{ii},xselected(jj,:));
            end
            Data.gevaltime=[Data.gevaltime;toc(gevalst)];
            Data.g=[Data.g;Gselected]; %update matrix with constraint function values
            Data.pen=[Data.pen; sum(max(zeros(1,length(Data.constraint)),Gselected).^2)];
        end
    end
    
    %check if new points are feasible and improvements
    for ii = 1:no_xnew
        if (Fselected(ii) < Data.fbest) %point has better function value then best point found so far
            if isfield(Data,'constraint') && Data.pen(end-no_xnew+ii)<=0 %better feasible point found
                Data.xbest = xselected(ii,:); %update best point found so far
                Data.fbest = Fselected(ii); %update best function value found so far
            elseif ~isfield(Data,'constraint') %no constraints
                Data.xbest = xselected(ii,:); %update best point found so far
                Data.fbest = Fselected(ii); %update best function value found so far
            end
        end
    end
    
    Data.S = [Data.S; xselected]; %update matrix of sample sites
    Data.Y = [Data.Y; Fselected]; %update vector of objective function values
    if isfield(Data,'constraint') %check which points feasible, and which are not
        infeas=find(Data.pen>0); %find indices of infeasible points
        feas=find(Data.pen <=0); %find indices of feasible points
    end
        
    if isfield(Data,'constraint') %if there are constraints, equip objective function value of infeasible points with penalty
        Yfeas=Data.Y(feas); %function values of feasible points
        Data.Ypenalized=zeros(length(Data.Y),1); 
        Data.Ypenalized(feas)=Data.Y(feas); %objective function values equipped with penalty
        Data.Ypenalized(infeas)= max(Yfeas)+100*Data.pen(infeas); %objective function values equipped with penalty
    else
        Data.Ypenalized=Data.Y; %if no constraints, objective function values stay as they are
    end
    
    Data.Ymed=Data.Ypenalized; %vector where large function values are replaced by median
    MedY=median(Data.Ypenalized); %median of function values 
    Data.Ymed(Data.Ypenalized>MedY)=MedY; %replace large function values by median, used in calculation of surrogate model parameters

    %update model parameters
    [lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate);

    save('Results','Data');
    fprintf('Number of function evaluation: %4.0f; Best feasible function value: %f\n', size(Data.S,1),Data.fbest)
end
end %function 