function CandValue=PredictFunctionValues(Data,Surrogate,CandPoint,lambda,gamma,...
    dmodel,beta,mmodel,w_m)

%predicts the objective function values using the desired surrogate model
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
%input
%Data - structure containing all information about the optimization problem
%Surrogate  - string with name of desired surrogate model
%CandPoint - matrix with candidates for next sample point. for these points
%the objective function values will be predicted
%lambda,gamma - parameters of RBF model, lambda=gamma=[] if RBF model not
%used
%dmodel - structure with parameters of kriging model, dmodel = [] if
%kriging model not used
%beta - parameters of polyniomial regression model, beta = [] if
%polynomial regression model not used
%mmodel - structure with parameters of MARS model, mmodel = [] if MARS
%model not used
%w_m - weights for models in mixtures, w_m = [] if mixture model not used
%
%Output
%CandValue - predicted objective function values for points in CandPoint
%--------------------------------------------------------------------------


%----------RBF models---------
if strcmp(Surrogate,'RBFlin') %linear RBF model
    CandValue=RBF_eval(CandPoint,Data.S,lambda,gamma,'linear'); 
elseif strcmp(Surrogate,'RBFtps') %thin plate spline RBF model
    CandValue=RBF_eval(CandPoint,Data.S,lambda,gamma,'TPS'); 
elseif strcmp(Surrogate,'RBFcub')  % cubic RBF model  
    CandValue=RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic'); 
%----------Kriging models---------
elseif strcmp(Surrogate(1:4),'KRIG') %all Kriging models, specific settings contained in structure dmodel
    CandValue = predictor(CandPoint, dmodel);        
%----------Polynomial models---------
elseif strcmp(Surrogate,'POLYlin')  %linear polynomial
    CandValue=POLY_eval(CandPoint,beta,'lin');
elseif strcmp(Surrogate,'POLYquad') %quadratic polynomial        
    CandValue=POLY_eval(CandPoint,beta,'quad'); 
elseif strcmp(Surrogate,'POLYquadr') %reduced quadratic polynomial        
    CandValue=POLY_eval(CandPoint,beta,'quadr');
elseif strcmp(Surrogate,'POLYcub')  %cubic polynomial               
    CandValue=POLY_eval(CandPoint,beta,'cub'); %not tested
elseif strcmp(Surrogate,'POLYcubr')  %reduced cubic polynomial                     
    CandValue=POLY_eval(CandPoint,beta,'cubr');
%----------MARS model---------    
elseif strcmp(Surrogate,'MARS') %Multivariate adaptive regression spline
    CandValue=arespredict(mmodel,CandPoint);
%----------Mixture models (2)---------    
elseif strcmp(Surrogate,'MIX_RcKg')  %Mixture model: cubic RBF & Kriging 
    CandValue_Rc=RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model for all candidate points
    CandValue_Kg = predictor(CandPoint, dmodel);        %objective function value prediction with kriging model for all candidate points
    CandValue= w_m(1)*CandValue_Rc+w_m(2)*CandValue_Kg; %weighted objective function value prediction
elseif strcmp(Surrogate,'MIX_RcM')  %Mixture model: cubic RBF & MARS
    CandValue_Rc=RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model for all candidate points
    CandValue_M=arespredict(mmodel,CandPoint); % objective function value prediction with MARS model for all candidate points
    CandValue= w_m(1)*CandValue_Rc+w_m(2)*CandValue_M; %weighted objective function value prediction
elseif strcmp(Surrogate,'MIX_RcPc')  %Mixture model: cubic RBF & MARS
    CandValue_Rc=RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic'); % objective function value prediction with cubic RBF model for all candidate points
    CandValue_Pcr=POLY_eval(CandPoint,beta,'cub');%objective function value prediction with reduced cubic polynomial for all candidate points
    CandValue= w_m(1)*CandValue_Rc+w_m(2)*CandValue_Pcr; %weighted objective function value prediction            
elseif strcmp(Surrogate,'MIX_KgM') %Mixture model: Kriging & MARS
    CandValue_Kg = predictor(CandPoint, dmodel);   %objective function value prediction with kriging model for all candidate points
    CandValue_M=arespredict(mmodel,CandPoint);    % objective function value prediction with MARS model for all candidate points
    CandValue= w_m(1)*CandValue_Kg+w_m(2)*CandValue_M; %weighted objective function value prediction
elseif strcmp(Surrogate,'MIX_KgPcr')   %Mixture model: Kriging & reduced cubic polynomial 
    CandValue_Kg = predictor(CandPoint, dmodel);  %objective function value prediction with kriging model for all candidate points
    CandValue_Pcr=POLY_eval(CandPoint,beta,'cubr');%objective function value prediction with reduced cubic polynomial for all candidate points
    CandValue= w_m(1)*CandValue_Kg+w_m(2)*CandValue_Pcr; %weighted objective function value prediction
elseif strcmp(Surrogate,'MIX_KgPc')   %Mixture model: Kriging & reduced cubic polynomial 
    CandValue_Kg = predictor(CandPoint, dmodel);  %objective function value prediction with kriging model for all candidate points
    CandValue_Pcr=POLY_eval(CandPoint,beta,'cub');%objective function value prediction with reduced cubic polynomial for all candidate points
    CandValue= w_m(1)*CandValue_Kg+w_m(2)*CandValue_Pcr; %weighted objective function value prediction
elseif strcmp(Surrogate,'MIX_KgPqr')  %Mixture model: Kriging & reduced quadratic polynomial   
    CandValue_Kg = predictor(CandPoint, dmodel);  %objective function value prediction with kriging model for all candidate points
    CandValue_Pqr=POLY_eval(CandPoint,beta,'quadr'); %objective function value prediction with reduced quadratic polynomial for all candidate points
    CandValue= w_m(1)*CandValue_Kg+w_m(2)*CandValue_Pqr; %weighted objective function value prediction
%add more surrogate models here if desired

%----------Mixture models (3)---------        
elseif strcmp(Surrogate,'MIX_RcKgM')   %Mixture model: cubic RBF &Kriging & reduced quadratic polynomial         
    CandValue_Rc=RBF_eval(CandPoint,Data.S,lambda,gamma,'cubic');  % objective function value prediction with cubic RBF model for all candidate points
    CandValue_Kg = predictor(CandPoint, dmodel);  %objective function value prediction with kriging model for all candidate points
    CandValue_M=arespredict(mmodel,CandPoint);  % objective function value prediction with MARS model for all candidate points
    CandValue= w_m(1)*CandValue_Rc+w_m(2)*CandValue_Kg + w_m(3)*CandValue_M; %weighted objective function value prediction
end

end%function