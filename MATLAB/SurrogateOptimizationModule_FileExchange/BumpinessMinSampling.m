function Data = BumpinessMinSampling(Data, maxeval, Surrogate, lambda, gamma, tolerance)

%samples points where bumpiness measure is minimized
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
%
%Input
%Data - structure, contains all information about the optimization problem
%maxeval - integer, maximum number of allowed function evaluations
%Surrogate - string, surrogate model type to be used
%lambda,gamma  - vectors, parameters of RBF model; lambda=gamma=[] if RBF model not used
%tolerance - scalar, distance when two points are considered equal
%
%Output
%Data - structure, contains updated information about the optimization problem
%--------------------------------------------------------------------------   


%different RBF models
if strcmp(Surrogate,'RBFlin')
    flag='linear'; %linear RBF
elseif strcmp(Surrogate,'RBFtps')
    flag='TPS';    %thin plate spline RBF
elseif strcmp(Surrogate,'RBFcub')  
    flag='cubic';%cubic RBF
end

iterctr=1; %initialize iteration counter, needed for computing objective function value target
w_j=[0 0.0001 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 ...
    0.11 0.12 0.13 0.15 0.20 0.25 0.30 0.40 0.50 0.75 1.0 1.5 2 3 100]; %weight factors used in global grid search,
%uses in every iteration a different weight for creating different
%target values; 
N=length(w_j); %cycle length, +1 because (N+1)th cycle is local search
w_counter=0;

while size(Data.S,1) < maxeval %iterate until stopping criterion met
    %first find minimum of response surface
    myfun=@(x)RBF_eval(x,Data.S,lambda,gamma,flag); % RBF model for predicting objective function values
    %use ordinary DDS algorithm for finding the minimum of the response surface
    %start point for DDS is randomly generated from variable domain
    x0=Data.xlow+rand(1,Data.dim).*(Data.xup-Data.xlow);
    [xmin , ymin]=ODDS(Data,x0,myfun, tolerance);


    %%%%%%%%%%%%%%%%%%%%------------------
    % add here other methods for finding the minimum of the response
    % measure
    %%%%%%%%%%%%%%%%%%%%--------------------


    %compute fdelta for computing target of objective function value
    if Data.fbest > 0
        fdelta=min( max( 1, Data.fbest ), max(Data.Ymed)-Data.fbest );
    else
        fdelta= min( 10*max(  1,abs(Data.fbest) ), max(Data.Ymed)-Data.fbest );  
    end

    if mod(iterctr,N)==0 %local search
        if ymin < Data.fbest-10e-6*abs(Data.fbest) %uses min of response surface as new sample site
            xselected=xmin;
            fevalst=tic; %do expensive objective function evaluation and record time needed
            Fselected = feval(Data.objfunction,xselected);
            Data.fevaltime=[Data.fevaltime;toc(fevalst)]; %record time for function evaluatin

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
            MedY=median(Data.Y);
            Data.Ymed(Data.Y>MedY)=MedY;                                 
            %and update the rsponse surface parameters
            [lambda, gamma]=RBF(Data.S,Data.Ymed,flag); 
            iterctr=iterctr+1;%update iteration number counter (needed for adjusting w_j)
            continue; %go to next iteration
        else
            ObjFtargets=Data.fbest-10e-2*abs(Data.fbest);  %determines target for objective function value for local search
        end          
    else 
        w_counter=w_counter+1;
        if mod(w_counter, N)==0
            w=w_j(end);
        else
            w=w_j(mod(w_counter, N)); %determine w_j for current objective function value target
        end
        ObjFtargets=ymin-w*fdelta;  %current objective function value target    
    end

    [lambda, gamma]=RBF(Data.S,Data.Ymed,flag); % compute model parameters         
    %define bumpiness measure function
    hn=@(x)bumpiness_measure(x,Data,flag, ObjFtargets, tolerance, lambda, gamma);
    %use ordinary DDS algorithm for finding the minimum of bumpiness
    %measure
    x0=Data.xlow+rand(1,Data.dim).*(Data.xup-Data.xlow);
    xselected=ODDS(Data,  x0,hn, tolerance);
    %%%%%%%%%%%%%%%%%%%%------------------
    % add here other methods for finding the minimum of the bumpiness
    % measure
    %%%%%%%%%%%%%%%%%%%%--------------------

    % perform function evaluation at the selected point
    fevalst=tic;
    Fselected = feval(Data.objfunction,xselected);
    Data.fevaltime=[Data.fevaltime;toc(fevalst)]; %record objective function evaluation time

    if (Fselected < Data.fbest) %new point is better than best point found so far?
        Data.xbest = xselected; %update best point found so far
        Data.fbest = Fselected;
    end

    %update data vectors
    Data.S = [Data.S; xselected]; %sample site matrix
    Data.Y = [Data.Y; Fselected]; %objective function values
    Data.Ymed = Data.Y; %data where large function values set to median, for calculation of surrogate model parameters
    % replace large function values by the median of all available function
    % values
    MedY=median(Data.Y);
    Data.Ymed(Data.Y>MedY)=MedY; 
    %and update the rsponse surface parameters
    [lambda, gamma]=RBF(Data.S,Data.Ymed,flag); 
    iterctr=iterctr+1; %update iteration counter (needed for adjusting w_j)
    fprintf('Number of function evaluation: %4.0f; Best feasible function value: %f\n', size(Data.S,1),Data.fbest)
    save('Results','Data');
end

end%function