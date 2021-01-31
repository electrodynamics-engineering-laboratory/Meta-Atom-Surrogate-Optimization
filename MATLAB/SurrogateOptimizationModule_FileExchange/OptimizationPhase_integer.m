function Data= OptimizationPhase_integer(Data,surrogate_model,maxeval)

%for problems where all variables have integer constraints. Problems may
%have additional constraints whose handles are in the fields
%Data.constraint{1}, Data.constraint{2}, ...
%If there are additional constraints and the initial experimental design
%does not contain any feasible points, optimizaion phase 1 tries to find a
%feasible point by minimizing a constraint violation function.
%After a feasible point has been found (or the maximal number of function
%evaluations has been reached), optimization phase 1 stops, and
%optimization phase 2 starts. The goal of optimization phase 2 is to find
%further feasible improvements of the best objective function value found
%so far. 
%Both optimization phases use the candidate point sampling strategy.
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
%Data - structure, contains all information about the optimization problem
%surrogate_model - string, contains the name of the desired (mixture)
%surogate model
%maxeval - maximum number of allowed function evaluations
%
%Output:
%Data - updated structure with all information about the problem
%--------------------------------------------------------------------------

%determine perturbation probability for each variable when generating
%candidate points
if Data.dim>5 %if problem dimension is larger than 5, perturb each variable with prob. 5/dimension
    P=max(0.1,5/Data.dim);
else %if problem dimension is <= 5, perturb all variables with prob. 1
    P=1;
end

stdev_int=[1,2,3]; %perturbation ranges 

%compute constraint function values at points in starting design
if isfield(Data,'constraint') %check if problem is constrained
    pen=zeros(size(Data.S,1),length(Data.constraint));
    Data.gevaltime=zeros(size(Data.Y));
    for jj=1:size(Data.S,1) %check for every point if constraints satisfied
        gevalst=tic; %record the time of constraint evaluations
        for ii = 1:length(Data.constraint)
            pen(jj,ii)=feval(Data.constraint{ii},Data.S(jj,:));
        end
        Data.gevaltime(ii,1)=toc(gevalst);
    end
    Data.g=pen; %matrix containing constraint function values
    % sum up squares of constraint function values that are larger than 0 
    Data.pen=sum(max(zeros(size(pen)),pen).^2,2); 
    if ~any(Data.pen <= 0) %no feasible point in initial experimental design
        Data.Feasible=false; 
        %optimization phase 1: minimize constraint violation function
        Data = SOI_OP1(Data,maxeval,P,stdev_int,surrogate_model);
    else
        Data.Feasible=true; %skip phase 1
        S_feas=Data.S(Data.pen<=0,:); %feasible points
        Y_feas=Data.Y(Data.pen<=0); %feasible function values
        [Data.fbest,idx]=min(Y_feas); %best feasible function value
        Data.xbest=S_feas(idx,:); %best feasible point
    end

    if ~Data.Feasible %if no feasible point found and max number function evaluations reached
        Data.xbest=NaN*ones(1,Data.dim);
        Data.fbest=NaN;
    end
else %there are no constraints
    [Data.fbest,idx]=min(Data.Y);
    Data.xbest=Data.S(idx,:);    
end
    
%optimization ohase 2 (if still function evaluations left)
if size(Data.S,1) < maxeval
    Data = SOI_OP2(Data,maxeval,P,stdev_int,surrogate_model);
end
end %function