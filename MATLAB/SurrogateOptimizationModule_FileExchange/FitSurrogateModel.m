function [lambda,gamma,dmodel,mmodel,beta, w_m] = FitSurrogateModel(Data, Surrogate)

%compute parameters of desired surrogate model(s)
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
%Surrogate - string containing surrogate mdoel to be used
%
%Output:
%lambda,gamma  - vectors, parameters of RBF model; lambda=gamma=[] if RBF model not used
%dmodel - structure, parameters for kriging model; dmodel=[] if kriging model not used
%mmodel - structure, parameters for MARS model; mmodel=[] if MARS model not used
%beta - vector, parameters of polynomial regression model; beta=[] if polynomial model not used
%w_m - vector, contains the weights for the models in mixtures; w_m=[] if no mixture model used
%--------------------------------------------------------------------------


%initialize parameters of surrogate models
lambda=[];
gamma=[];
dmodel=[];
mmodel=[];
beta=[];
w_m=[];

%---------RBF models----------
if strcmp(Surrogate,'RBFlin') %linear RBF model
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'linear'); 
elseif strcmp(Surrogate,'RBFtps') %thin plate spline RBF model
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'TPS'); 
elseif strcmp(Surrogate,'RBFcub')    %cubic RBF model
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); 
%---------Kriging models uses DACE----------    
elseif strcmp(Surrogate,'KRIGcub0')  %kriging model with cubic correlation function and 0th order regression polynomial  
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly0', 'corrcubic', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
elseif strcmp(Surrogate,'KRIGcub1')  %kriging model with cubic correlation function and 1st order regression polynomial    
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'corrcubic', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
elseif strcmp(Surrogate,'KRIGcub2')  %kriging model with cubic correlation function and 2nd order regression polynomial    
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly2', 'corrcubic', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
elseif strcmp(Surrogate,'KRIGexp0')  %kriging model with exponential correlation function and 0th order regression polynomial    
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly0', 'correxp', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
elseif strcmp(Surrogate,'KRIGexp1')  %kriging model with exponential correlation function and 1st order regression polynomial      
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'correxp', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
elseif strcmp(Surrogate,'KRIGexp2')  %kriging model with exponential correlation function and 2nd order regression polynomial      
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly2', 'correxp', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));      
elseif strcmp(Surrogate,'KRIGgexp0') %kriging model with generalized exponential correlation function and 0th order regression polynomial       
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly0', 'correxpg', ones(1,Data.dim+1), 1e-1*ones(1,Data.dim+1), [20*ones(1,Data.dim),2]);
elseif strcmp(Surrogate,'KRIGgexp1') %kriging model with generalized exponential correlation function and 1st order regression polynomial          
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'correxpg', ones(1,Data.dim+1), 1e-1*ones(1,Data.dim+1), [20*ones(1,Data.dim),2]);  
elseif strcmp(Surrogate,'KRIGgexp2') %kriging model with generalized exponential correlation function and 2nd order regression polynomial          
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly2', 'correxpg', ones(1,Data.dim+1), 1e-1*ones(1,Data.dim+1), [20*ones(1,Data.dim),2]);     
elseif strcmp(Surrogate,'KRIGgauss0') %kriging model with gaussian correlation function and 0th order regression polynomial       
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly0', 'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
elseif strcmp(Surrogate,'KRIGgauss1') %kriging model with gaussian correlation function and 1st order regression polynomial          
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
elseif strcmp(Surrogate,'KRIGgauss2') %kriging model with gaussian correlation function and 2nd order regression polynomial          
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly2', 'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
elseif strcmp(Surrogate,'KRIGlin0')   %kriging model with linear correlation function and 0th order regression polynomial            
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly0', 'corrlin', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
elseif strcmp(Surrogate,'KRIGlin1')   %kriging model with linear correlation function and 1st order regression polynomial             
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'corrlin', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
elseif strcmp(Surrogate,'KRIGlin2')   %kriging model with linear correlation function and 2nd order regression polynomial             
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly2', 'corrlin', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
elseif strcmp(Surrogate,'KRIGspline0') %kriging model with spline correlation function and 0th order regression polynomial                   
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly0', 'corrspline', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
elseif strcmp(Surrogate,'KRIGspline1') %kriging model with spline correlation function and 1st order regression polynomial                      
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'corrspline', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
elseif strcmp(Surrogate,'KRIGspline2') %kriging model with spline correlation function and 2nd order regression polynomial                      
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly2', 'corrspline', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
elseif strcmp(Surrogate,'KRIGsphere0')  %kriging model with spherical correlation function and 0th order regression polynomial                   
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly0', 'corrspherical', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
elseif strcmp(Surrogate,'KRIGsphere1')  %kriging model with spherical correlation function and 1st order regression polynomial                     
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'corrspherical', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
elseif strcmp(Surrogate,'KRIGsphere2')  %kriging model with spherical correlation function and 2nd order regression polynomial                     
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly2', 'corrspherical', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
%---------polynomial regression models----------    
elseif strcmp(Surrogate,'POLYlin')   %linear polynomial model 
    beta=POLY(Data.S,Data.Ymed,'lin');
elseif strcmp(Surrogate,'POLYquad')  %quadratic polynomial model        
    beta=POLY(Data.S,Data.Ymed,'quad'); 
elseif strcmp(Surrogate,'POLYquadr') %reduced quadratic polynomial model        
    beta=POLY(Data.S,Data.Ymed,'quadr');
elseif strcmp(Surrogate,'POLYcub')   %cubic polynomial model        
    beta=POLY(Data.S,Data.Ymed,'cub');
elseif strcmp(Surrogate,'POLYcubr')  %reduced cubic polynomial model               
    beta=POLY(Data.S,Data.Ymed,'cubr');
%---------MARS model----------    
elseif strcmp(Surrogate,'MARS') %multivariate adaptive regression spline (MARS)
    mmodel = aresbuild(Data.S,Data.Ymed);
    
%------------------------------
% add more models here if necessary. Alter accordingly later parts in the
% code
%------------------------------

%---------Mixture models (2)----------    
elseif strcmp(Surrogate,'MIX_RcKg') %mixture model of cubic RBF and Kriging model with gaussian correlation and 1st order regression polynomial
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); % RBF
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  %Kriging
    %w_m=[.5, .5]; %set fixed weights if desired (sum must be one)
    w_m=DempsterFor2models(Data,Surrogate); %uses dempster shafer theory to adjust model weights
elseif strcmp(Surrogate,'MIX_RcM') %mixture model of cubic RBF and MARS
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); % RBF
    mmodel = aresbuild(Data.S,Data.Ymed); %MARS
    %w_m=[.5, .5]; %set fixed weights if desired (sum must be one)
    w_m=DempsterFor2models(Data,Surrogate); %uses dempster shafer theory to adjust model weights
elseif strcmp(Surrogate,'MIX_RcPc') %mixture model of Kriging model with gaussian correlation and 1st order regression polynomial and reduced cubic polynomial
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); % RBF
    beta=POLY(Data.S,Data.Ymed,'cub'); %Polynomial
    %w_m=[.5, .5]; %set fixed weights if desired (sum must be one)
    w_m=DempsterFor2models(Data,Surrogate); %uses dempster shafer theory to adjust model weights
elseif strcmp(Surrogate,'MIX_KgM') %mixture model of Kriging model with gaussian correlation and 1st order regression polynomial and MARS
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  %Kriging
    mmodel = aresbuild(Data.S,Data.Ymed); %MARS
    %w_m=[.5, .5]; %set fixed weights if desired (sum must be one)
    w_m=DempsterFor2models(Data, Surrogate); %uses dempster shafer theory to adjust model weights
elseif strcmp(Surrogate,'MIX_KgPcr') %mixture model of Kriging model with gaussian correlation and 1st order regression polynomial and reduced cubic polynomial
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  %Kriging
    beta=POLY(Data.S,Data.Ymed,'cubr'); %reduced cubic polynomial
    %w_m=[.5, .5]; %set fixed weights if desired (sum must be one)
    w_m=DempsterFor2models(Data, Surrogate); %uses dempster shafer theory to adjust model weights
elseif strcmp(Surrogate,'MIX_KgPc') %mixture model of Kriging model with gaussian correlation and 1st order regression polynomial and reduced cubic polynomial
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  %Kriging
    beta=POLY(Data.S,Data.Ymed,'cub'); %Polynomial
    %w_m=[.5, .5]; %set fixed weights if desired (sum must be one)
    w_m=DempsterFor2models(Data, Surrogate); %uses dempster shafer theory to adjust model weights
elseif strcmp(Surrogate,'MIX_KgPqr') %mixture model of Kriging model with gaussian correlation and 1st order regression polynomial and reduced quadratic polynomial
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  %Kriging
    beta=POLY(Data.S,Data.Ymed,'quadr'); %reduced quadratic polynomial
    %w_m=[.5, .5]; %set fixed weights if desired (sum must be one)
    w_m=DempsterFor2models(Data, Surrogate); %uses dempster shafer theory to adjust model weights
    
% add more 2-model combinations here if necessary. Update function "DempsterFor2models" accordingly 

%---------Mixture models (3)----------        
elseif strcmp(Surrogate,'MIX_RcKgM') %mixture model of cubic RBF, Kriging model with gaussian correlation and 1st order regression polynomial and MARS
    [lambda, gamma]=RBF(Data.S,Data.Ymed,'cubic'); %  RBF
    dmodel = dacefit(Data.S, Data.Ymed, 'regpoly1', 'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));   %Kriging
    mmodel = aresbuild(Data.S,Data.Ymed); %MARS
    %w_m=[1/3, 1/3, 1/3]; %set fixed weights if desired (sum must be one)
    w_m=DempsterFor3models(Data,Surrogate); %uses dempster shafer theory to adjust model weights
    
% add more 3-model combinations here if necessary. Update function "DempsterFor3models" accordingly     
    
end

end%function