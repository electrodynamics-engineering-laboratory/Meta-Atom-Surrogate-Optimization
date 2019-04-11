function Data = SurfaceMinSampling(Data,maxeval,Surrogate,lambda,gamma,...
        dmodel,mmodel,beta, w_m, tolerance)
    
%uses minimum point of response surface to determinte the next sample point in each iteration    
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

while size(Data.S,1) < maxeval  
    %determine which respose surface to minimize
    %------------RBF models ----------------
    if strcmp(Surrogate,'RBFlin') %linear RBF model
        myfun=@(x)RBF_eval(x,Data.S,lambda,gamma,'linear'); 
    elseif strcmp(Surrogate,'RBFtps') %thin plate spline RBF model
        myfun=@(x)RBF_eval(x,Data.S,lambda,gamma,'TPS'); 
    elseif strcmp(Surrogate,'RBFcub')  % cubic RBF model  
        myfun=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); 
    %------------Kriging models --------------    
    elseif strcmp(Surrogate(1:4),'KRIG')  %all Kriging models
        myfun=@(x)predictor(x, dmodel);        
    %------------Polynomial models --------------     
    elseif strcmp(Surrogate,'POLYlin')  %linear polynomial model
        myfun=@(x)POLY_eval(x,beta,'lin');
    elseif strcmp(Surrogate,'POLYquad')  %quadratic polynomial model       
        myfun=@(x)POLY_eval(x,beta,'quad'); 
    elseif strcmp(Surrogate,'POLYquadr') %reduced quadratic polynomial model
        myfun=@(x)POLY_eval(x,beta,'quadr');
    elseif strcmp(Surrogate,'POLYcub')  %cubic polynomial model       
        myfun=@(x)POLY_eval(x,beta,'cub'); 
    elseif strcmp(Surrogate,'POLYcubr')  %reduced cubic polynomial model      
        myfun=@(x)POLY_eval(x,beta,'cubr');
    %------------MARS models --------------     
    elseif strcmp(Surrogate,'MARS') %multivariate adaptive regression splines
        myfun=@(x)arespredict(mmodel,x);
    %------------Mixture models (2)--------------     
    elseif strcmp(Surrogate,'MIX_RcKg')   %mixture model of cubic RBF and Kriging with gaussian correlation and 1st order regression polynomial
        myfun1=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); % cubic RBF model
        myfun2=@(x)predictor(x, dmodel);  %kriging model
        myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted mixture
    elseif strcmp(Surrogate,'MIX_RcM') %mixture model of cubic RBF and MARS
        myfun1=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); % cubic RBF model
        myfun2=@(x)arespredict(mmodel,x); %MARS model
        myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted mixture
    elseif strcmp(Surrogate,'MIX_RcPc')  %Mixture model: cubic RBF & MARS
        myfun1=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model for all candidate points
        myfun2=@(x)POLY_eval(x,beta,'cub');%objective function value prediction with reduced cubic polynomial for all candidate points
        myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted objective function value prediction            
    elseif strcmp(Surrogate,'MIX_KgM') %mixture model of Kriging with gaussian correlation and 1st order regression polynomial and MARS model
        myfun1=@(x)predictor(x, dmodel);   %kriging model
        myfun2=@(x)arespredict(mmodel,x);  %MARS model
        myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted mixture
    elseif strcmp(Surrogate,'MIX_KgPcr') %mixture model of Kriging with gaussian correlation and 1st order regression polynomial, and reduced cubic polynomial
        myfun1=@(x)predictor(x, dmodel);     %kriging model
        myfun2=@(x)POLY_eval(x,beta,'cubr'); %reduced cubic model
        myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted mixture
    elseif strcmp(Surrogate,'MIX_KgPc') %mixture model of Kriging with gaussian correlation and 1st order regression polynomial, and  cubic polynomial
        myfun1=@(x)predictor(x, dmodel);     %kriging model
        myfun2=@(x)POLY_eval(x,beta,'cub'); % cubic model
        myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted mixture
    elseif strcmp(Surrogate,'MIX_KgPqr') %mixture model of Kriging with gaussian correlation and 1st order regression polynomial, and reduced quadratic polynomial   
        myfun1=@(x)predictor(x, dmodel);        %kriging model
        myfun2=@(x)POLY_eval(x,beta,'quadr');   %reduced quadratic model
        myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x); %weighted mixture
    %add other 2-model mixtures here    
    %------------Mixture models (3)--------------         
    elseif strcmp(Surrogate,'MIX_RcKgM') %mixture model of Kriging with gaussian correlation and 1st order regression polynomial, and reduced cubic polynomial       
        myfun1=@(x)RBF_eval(x,Data.S,lambda,gamma,'cubic'); % cubic RBF model
        myfun2=@(x)predictor(x, dmodel);        %kriging model
        myfun3=@(x)arespredict(mmodel,x); %MARS model
        myfun=@(x)w_m(1)*myfun1(x)+w_m(2)*myfun2(x)+ w_m(3)*myfun3(x); %weighted mixture
    end



    %option_1: ordinary DDS approach for minimization
    x0=Data.xlow+rand(1,Data.dim).*(Data.xup-Data.xlow); %initial  point for DDS
    xselected =ODDS(Data, x0,myfun, tolerance);


    %%%%%%%%%%%%%%%%%%%%%%%%----------------------
    % include different optimization routines for finding the minimum
    % of the response surface, eg fmincon
    %%%%%%%%%%%%%%%%%%%%%%%%----------------------


    % perform function evaluation at the selected point
    fevalst=tic;
    Fselected = feval(Data.objfunction,xselected);
    Data.fevaltime=[Data.fevaltime;toc(fevalst)]; %record time needed for function evaluation

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
    %recompute model parameters    
    [lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate);
    fprintf('Number of function evaluation: %4.0f; Best feasible function value: %f\n', size(Data.S,1),Data.fbest)
    save('Results','Data');
end
end%function