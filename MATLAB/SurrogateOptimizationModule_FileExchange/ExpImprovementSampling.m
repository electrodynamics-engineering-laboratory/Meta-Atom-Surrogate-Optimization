function Data = ExpImprovementSampling(Data,maxeval) 

%uses point that maximizes the expected improvement as new sample site
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
%Data - structure containing all information about the problem
%maxeval - maximum number of allowed function evaluations
%
%Output:
%Data - structure containing updated information about the problem
%--------------------------------------------------------------------------

while size(Data.S,1) < maxeval %do until stopping criterion met
    Data.UpperTheta=ones(1,Data.dim).*2; %upper bounds on model parameter theta
    Data.LowerTheta=ones(1,Data.dim).*-3; %lower bounds on model parameter theta
    [Data.dist,Data.ij]=distanceupdate(Data); %compute pairwise distances between all points
    optionsEI=gaoptimset('PopulationSize',10,'Display','off'); %set options for genetic algorithm (GA) that computes model parameters theta
    [Data.Theta,~] = ga(@(x)likelihood_new(Data,x),Data.dim,[],[],[],[],Data.LowerTheta,Data.UpperTheta,[],optionsEI); %optimizaes the theta's with GA
    [~,Data.Psi,Data.U]=likelihood_new(Data,Data.Theta); %uses optimized theta's and computes the parameters needed for computing expected improvement
    myEI=@(x) -log10(ExpectedImprovement(Data,x)) ; %negative log expected improvement function to be minimized
    %use the GA approach (used by forrester et al) to optimize expected
    %improvement
    xselected=ga(myEI,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],optionsEI);

    %%%%%%%%%%%%%%%%%%%%%%%%%--------------------
    % add here other options for finding the point that maximizes the
    % expected improvement
    %%%%%%%%%%%%%%%%%%%%%%%%%--------------------

    % perform function evaluation at the selected point
    fevalst=tic;
    Fselected = feval(Data.objfunction,xselected);
    Data.fevaltime=[Data.fevaltime;toc(fevalst)];%record time needed for the function evaluation
    if (Fselected < Data.fbest) %new point is better than best point found so far?
        Data.xbest = xselected; %update best point found so far
        Data.fbest = Fselected; %update best function value found so far
    end  
    %update data vectors
    Data.S = [Data.S; xselected]; %sample site matrix
    Data.Y = [Data.Y; Fselected]; %objective function values
    Data.Ymed = Data.Y; %data where large function values set to median, for calculation of surrogate model parameters
    % replace large function values by the median of all available function
    % values
    %MedY=median(Data.Y);
    %Data.Ymed(Data.Y>MedY)=MedY;        
    fprintf('Number of function evaluation: %4.0f; Best feasible function value: %f\n', size(Data.S,1),Data.fbest)
    save('Results','Data');
end

end %function