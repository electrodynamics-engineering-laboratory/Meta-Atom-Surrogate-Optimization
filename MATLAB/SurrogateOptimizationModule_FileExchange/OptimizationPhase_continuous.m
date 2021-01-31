function Data= OptimizationPhase_continuous(Data,Surrogate,SampleStrategy,maxeval)

%optimization phase for continuous unconstrained problems: 
%iteratively pick a new sample site, do the expensive
%function evaluation, and update the chosen surrogate model; iterate until
%stopping criterion met
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
%input:
%Data: structure array containing the optimization problem information
%Surrogate: string determining which surrogate model is being used
%SampleStrategy: string determining which strategy for picking the next
%                sample site is being used
%maxeval: maximum number of function evaluations allowed
%
%output: updated structure array "Data" containing all information about the
%solution
%--------------------------------------------------------------------------

tolerance = 1e-6;  % tolerance for deciding when two points coincide

[Data.fbest,idx]=min(Data.Y);  %best function value so far
Data.xbest=Data.S(idx,:); %best point found so far

Data.Ymed=Data.Y; %vector where large function values are replaced by median
MedY=median(Data.Y); %median of function values 
Data.Ymed(Data.Y>MedY)=MedY; %replace large function values by median, used in calculation of surrogate model parameters

% sampling strategy
if strcmp(SampleStrategy,'CAND') %candidate point sampling strategy  
    % determine parameters of surrogate model
    [lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate);
    Data = CandidatePointSampling(Data,maxeval,Surrogate,lambda,gamma,...
        dmodel,mmodel,beta, w_m, tolerance); 
elseif strcmp(SampleStrategy,'SURFmin') %minimum point of the response surface 
    % determine parameters of surrogate model
    [lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate);
    Data = SurfaceMinSampling(Data,maxeval,Surrogate,lambda,gamma,...
        dmodel,mmodel,beta, w_m, tolerance);    
%use the point that maximizes the expected improvement (only for kriging, 
%generalized exponential) parts of the codes are from Forrester et al book 
%for EGO codes options for genetic algorithm    
elseif strcmp(SampleStrategy,'EImax')   
    Data = ExpImprovementSampling(Data,maxeval);
%uses the point that minimizes a bumpiness function given a target for the objective function value    
%only for radial basis functions
%switches between local and global search
%determine which radial basis function model should be used
elseif strcmp(SampleStrategy,'BUMPmin') 
    % determine parameters of surrogate model
    [lambda,gamma,~,~,~, ~] = FitSurrogateModel(Data, Surrogate);
    Data = BumpinessMinSampling(Data, maxeval, Surrogate, lambda, gamma, tolerance);
%minimize scoring function          
elseif strcmp(SampleStrategy,'SCOREmin')
    %uses the minimum site of a scoring function as next sample site; this
    %corresponds to choosing the point with "optimal" score, as opposed to
    %the candidate point approach where the best point amongst a limited
    %number of candidate points is chosen
    % determine parameters of surrogate model
    [lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate);
    Data = ScoreMinSampling(Data, maxeval, Surrogate, lambda,gamma, dmodel,...
        mmodel, beta,w_m,tolerance);
end %selection of sampling strategy

end %function