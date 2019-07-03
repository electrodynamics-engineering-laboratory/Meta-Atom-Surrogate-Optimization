function [ys,dmodel] = get_ys_krig(x,y,xs,nvar,b,c,cor,reg)
% This function builds an PRS metamodel and predicts the function value at
% test points.

% ---Variable Descriptions---
% x         = normalized design(training points): in an 
%           [num_of_points by num_of_variables] matrix
% y         = vector of function values at each x
% xs        = test points: in an 
%           [num_of_points by num_of_variables] matrix
% nvar      = the number of design variables
% b         = the lower bound on theta 0.1
% c         = the upper bound on theta 0.9
% cor       = the correlation functions
%       --Supported Correlation Functions--
%    cor = 1 -> Cubic;
%    cor = 2 -> Exponential;
%    cor = 3 -> Gaussian;
%    cor = 4 -> Linear;
%    cor = 5 -> Spherical;
%    cor = 6 -> Spline;
% reg       = regression polynomial
%       --Supported Regression Polynomials--
%    reg = 1 -> 0 degree;
%    reg = 2 -> 1st degree;
%    reg = 3 -> 2nd degree;


% Add the path name where the dace toolbox folder is located
addpath('...\dace')

% Initial estimate of theta0 (correlation parameter)
a = (b+c)/2;

if cor == 1
    correlation = @corrcubic;
elseif cor == 2
    correlation = @correxp;
elseif cor == 3
    correlation = @corrgauss;
elseif cor == 4
    correlation = @corrlin;
elseif cor == 5
    correlation = @corrspherical;
elseif cor == 6
    correlation = @corrspline;
end
if reg == 1
    regression = @regpoly0;
elseif reg == 2
    regression = @regpoly1;
elseif reg == 3
    regression = @regpoly2;
end
% MODEL IDENTIFICATION 
% Identifying the Kriging Model Parameters (Correlation Structure)
% Initial values for the correlation parameters and lower and upper bounds
theta0 = a*ones(1,nvar); 
lob    = b*ones(1,nvar); 
upb    = c*ones(1,nvar);

% Specification of the input and output data, and the regression 
% and correlation model 
[dmodel, perf] = dacefit(x, y, regression, correlation, theta0, lob, upb);
% dmodel: a data structure with all the necessary information to feed the
% predictor function
% perf: a data structure with details of the optimization process leading
%to the correlation structure parameters  

% PREDICTION at test points
[ys, mse] = predictor(xs, dmodel);

% ys:  values of the prediction at the xs locations
% mse: mse for each predicted point