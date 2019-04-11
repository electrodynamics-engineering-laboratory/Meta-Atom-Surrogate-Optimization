function SurrogateModelModule_v1(data_file,maxeval,surrogate_model, sampling_technique, ...
    initial_design, number_startpoints, starting_point)

%Surrogate model toolbox for unconstrained continuous, constrained integer
%and constrained mixed-integer global optimization problems. 
%The user can choose beween different options for
%   - the surrogate model
%   - the sampling strategy
%   - the initial experimental design
%The user can determine the maximum number of allowed function evaluations.
%the number of points in the initial starting design, and one or more
%points that are added to the starting design

%--------------------------------------------------------------------------
%Copyright (c) 2012 by Juliane Mueller
%
% This file is part of the surrogate model module toolbox.
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%Tampere University of Technology, Finland
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------

%INPUT - mandatory
%data_file: string with name of data file containing optimization problem
%data
%INPUT - optional
%maxeval - maximum number of allowed function evaluations (default 400)
%surrogate_model: string defining which surrogate model to be used:
%   RBFcub - cubic RBF (default)
%   RBFtps - thin plate spline RBF
%   RBFlin - linear RBF
%   KRIGexp0 - Kriging with exponential correlation, regression polynomial
%               order 0
%   KRIGexp1 - Kriging with exponential correlation, regression polynomial
%               order 1
%   KRIGexp2 - Kriging with exponential correlation, regression polynomial
%               order 2
%   KRIGgexp0 - Kriging with generalized exponential correlation, regression polynomial
%               order 0
%   KRIGgexp1 - Kriging with generalized exponential correlation, regression polynomial
%               order 1
%   KRIGgexp2 - Kriging with generalized exponential correlation, regression polynomial
%               order 2
%   KRIGgauss0 - Kriging with Gaussian correlation, regression polynomial
%               order 0
%   KRIGgauss1 - Kriging with Gaussian correlation, regression polynomial
%               order 1
%   KRIGgauss2 - Kriging with Gaussian correlation, regression polynomial
%               order 2
%   KRIGlin0 - Kriging with linear correlation, regression polynomial
%               order 0
%   KRIGlin1 - Kriging with linear correlation, regression polynomial
%               order 1
%   KRIGlin2 - Kriging with linear correlation, regression polynomial
%               order 2
%   KRIGspline0 - Kriging with spline correlation, regression polynomial
%               order 0
%   KRIGspline1 - Kriging with spline correlation, regression polynomial
%               order 1
%   KRIGspline2 - Kriging with spline correlation, regression polynomial
%               order 2
%   KRIGsphere0 - Kriging with spherical correlation, regression polynomial
%               order 0
%   KRIGsphere1 - Kriging with spherical correlation, regression polynomial
%               order 1
%   KRIGsphere2 - Kriging with spherical correlation, regression polynomial
%               order 2
%   KRIGcub0 - Kriging with cubic correlation, regression polynomial
%               order 0
%   KRIGcub1 - Kriging with cubic correlation, regression polynomial
%               order 1
%   KRIGcub2 - Kriging with cubic correlation, regression polynomial
%               order 2
%   POLYlin - linear polynomial
%   POLYquad - quadratic polynomial
%   POLYquadr - reduced quadratic polynomial
%   POLYcub - cubic polynomial
%   POLYcubr - reduced cubic polynomial
%   MARS - multivariate adaptive regression spline
%   MIX_RcKg - mixture model cubic RBF and gaussian Kriging
%   MIX_RcM - mixture model cubic RBF and MARS
%   MIX_RcPc - mixture model cubic RBF and cubic polynomial
%   MIX_KgM - mixture model gaussian Kriging and MARS
%   MIX_KgPc - mixture model gaussian Kriging and cubic polynomial
%   MIX_KgPqr - mixture model gaussian Kriging and reduced quadratic polynomial
%   MIX_KgPcr - mixture model gaussian Kriging and reduced cubic polynomial
%   MIX_RcKgM - mixture model cubic RBF, gaussian Kriging and MARS
%   ...define further mixtures if needed...

%sampling_technique: string, determining the technique for selecting the 
%next sample site
%   CAND - candidate poind approach (default)
%   SURFmin - uses the minimum point of response surface
%   EImaxmin - uses the point that maximizes the expected improvement
%       (currently working only for Kriging with generalized exponential
%       correlation function - implementation by Forrester)
%   SCOREmin - uses the scoring criteria from the candidate point approach
%       but tries to find the best point in the whole domain (instead of 
%       choosing the best point among a limited number of points as in
%       CAND)
%   BUMPmin - minimizes bumpiness measure (Gutmann, 2001), only for RBFs

%initial_design: string defining which initial experimental design should
%be used
%   SLHD - symmetric latin hypercube design 
%   LHS  - Matlab's latin hypercube design (default)
%   CORNER - corner points
%   SPACEFIL - space filling design (as in EGO)

%number_startpoints: integer defining the number of points in the initial
%experimental design, default 2(d+1)

%starting_point - matrix of dimension (m x d) where d is problem dimension,
%contains points that are added to the initial experimental design

%OUTPUT: none, results are saved in file Results.mat
%--------------------------------------------------------------------------
warning off
timeStart=tic; %record time of algorithm

if nargin < 1 ||  ~ischar(data_file) %check if input data file present
    error('You have to supply a file name with your data. See example files for information how to define problems.')
end
Data=feval(data_file); % load problem data
%check if all upper and lower bounds are given
if length(Data.xlow) ~= Data.dim || length(Data.xup) ~= Data.dim
    error('Vector length of lower and upper bounds must equal problem dimension')
end
%check if lower bounds < upper bounds
if any(Data.xlow >Data.xup)
    error('Lower bounds have to be lower than upper bounds.')
end

%check how many input arguments are given and use default values  if
%necessary
if nargin<2 %default settings for everything
    fprintf(['No surrogate model, sampling strategy, initial design, \n maximal number of allowed function evaluations, ',...
    'and number of \n starting points defined. Using default values.\n']);
    surrogate_model='RBFcub';   %if no surrogate model defined, use cubic RBF
    sampling_technique='CAND';  %if no sampling method definded, use candidate points
    initial_design='LHS';       %if no initial experimental design strategy defined, use symm. Latin hypercube
    maxeval= 400;               %if maximum number of function evaluations not defined, set to 400
    number_startpoints=2*(Data.dim+1); %default is number of starting points not defined
elseif nargin < 3
    fprintf(['No surrogate model, sampling strategy, initial design,\n and number of starting points '...
        'defined.\n Using default values.\n']);
    surrogate_model='RBFcub';   %if no surrogate model defined, use cubic RBF
    sampling_technique='CAND';  %if no sampling method definded, use candidate points
    initial_design='LHS';      %if no initial experimental design strategy defined, use symm. Latin hypercube    number_startpoints=2*(dim+1); %add the (2d+2)th point from list of starting points
    number_startpoints=2*(Data.dim+1); %default is number of starting points not defined
elseif nargin < 4
    fprintf(['No initial design strategy, sampling strategy,\n ',...
    'and number of starting points defined.\n Using default values.\n']);
    sampling_technique='CAND';  %if no sampling method definded, use candidate points
    initial_design='LHS';      %if no initial experimental design strategy defined, use symm. Latin hypercube    number_startpoints=2*(dim+1); %add the (2d+2)th point from list of starting points
    number_startpoints=2*(Data.dim+1); %default is number of starting points not defined
elseif nargin < 5
    fprintf(['No initial design strategy,\n  ',...
    'and number of starting points defined. Using default values.\n']);
    if strcmp(sampling_technique,'EImax') %special case for expected improvement sampling technique
        fprintf(['Using the expected improvement as criterion for selecting the next \n sample site works only with Kriging ' ...
            '- generalized exponential correlation;\n I will use this surrogate model!\n'])
       surrogate_model='KRIGgexp1';
    elseif strcmp(sampling_technique,'BUMPmin') %special case for bumpiness measure sampling technique
        if ~strcmp(surrogate_model,'RBFcub') || ~strcmp(surrogate_model,'RBFtps') ||...
                ~strcmp(surrogate_model,'RBFlin')
            fprintf(['Using the minimum of the bumpiness measure as criterion for selecting \n the next sample site requires using cubic, linear or thin plate spline  RBF\n' ...
            '- I will use cubic RBF from here!\n'])
            surrogate_model='RBFcub';   
        end
    end
    initial_design='LHS';
    number_startpoints=2*(Data.dim+1); %default is number of starting points not defined
elseif nargin < 6
    fprintf('No number of starting points defined. Using default value.\n');
    if strcmp(sampling_technique,'EImax')
        fprintf(['Using the expected improvement as criterion for picking the next sample site \n works only with Kriging \n' ...
            '- generalized exponential correlation; I will use this currogate model!\n'])
       surrogate_model='KRIGgexp1';
    elseif strcmp(sampling_technique,'BUMPmin') %special case for bumpiness measure sampling technique
        if ~strcmp(surrogate_model,'RBFcub') || ~strcmp(surrogate_model,'RBFtps') ||...
                ~strcmp(surrogate_model,'RBFlin')
            fprintf(['Using the minimum of the bumpiness measure as criterion for picking the next \n sample site requires using cubic, linear or thin plate spline  RBF\n' ...
            '- I will use cubic RBF from here!\n'])
            surrogate_model='RBFcub';   
        end
    end
    number_startpoints=2*(Data.dim+1); %default is number of starting points not defined
elseif nargin < 7
    if strcmp(sampling_technique,'EImax')
        fprintf(['Using the expected improvement as criterion for picking the next sample \n site works only with Kriging \n' ...
            '- generalized exponential correlation; I will use this currogate model!\n'])
        surrogate_model='KRIGgexp1';
    elseif strcmp(sampling_technique,'BUMPmin') %special case for bumpiness measure sampling technique
        if ~strcmp(surrogate_model,'RBFcub') || ~strcmp(surrogate_model,'RBFtps') ||...
                ~strcmp(surrogate_model,'RBFlin')
            fprintf(['Using the minimum of the bumpiness measure as criterion for picking the next sample \n site requires using cubic, linear or thin plate spline  RBF \n' ...
            '- I will use cubic RBF from here!\n'])
            surrogate_model='RBFcub';   
        end    
    end   
end

%check if starting point is of right dimension
if exist('starting_point','var') && ~isempty(starting_point)
    [m_sp,n_sp]=size(starting_point);
    if n_sp ~= Data.dim
        if m_sp == Data.dim
            starting_point=starting_point';
        else
            error('Provided starting point has incorrect dimension.')
        end
    end
end

%check if there are enough starting points for using RBF models
if any(strcmp(surrogate_model,{'RBFcub', 'RBFlin','RBFtsp'})) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) <Data.dim+1)
        fprintf('With RBF model I need at least Data.dim+1 starting points.\n')
        number_startpoints =Data.dim+1-size(starting_point,1);
    elseif(number_startpoints <Data.dim+1)  
        fprintf('With RBF model I need at least Data.dim+1 starting points.\n')
        number_startpoints =Data.dim+1;
    end
end

%check if there are enough starting points for using Kriging models with
%zero-order regression polynomials
if any(strcmp(surrogate_model,{'KRIGexp0', 'KRIGgexp0','KRIGgauss0','KRIGlin',...
        'KRIGspline0','KRIGsphere0','KRIGcub0'})) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) <Data.dim+1)
        fprintf('With Kriging and zeroth order polynomial I need at least Data.dim+1 starting points.\n')
        number_startpoints =Data.dim+1-size(starting_point,1);
    elseif(number_startpoints <Data.dim+1)  
        fprintf('With Kriging and zeroth order polynomial I need at least Data.dim+1 starting points.\n')
        number_startpoints =Data.dim+1;
    end
end

%check if there are enough starting points for using Kriging models with
%first-order regression polynomials
if any(strcmp(surrogate_model,{'KRIGexp1', 'KRIGgexp1','KRIGgauss1','KRIGlin1',...
        'KRIGspline1','KRIGsphere1','KRIGcub1'})) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) <Data.dim+1) 
        fprintf('With Kriging and first order polynomial I need at least Data.dim+1 starting points.\n')
        number_startpoints =Data.dim+1-size(starting_point,1);
    elseif(number_startpoints <Data.dim+1)  
        fprintf('With Kriging and first order polynomial I need at least Data.dim+1 starting points.\n')
        number_startpoints =Data.dim+1;
    end
end

%check if there are enough starting points for using Kriging models with
%second-order regression polynomials
if any(strcmp(surrogate_model,{'KRIGexp2', 'KRIGgexp2','KRIGgauss2','KRIGlin2',...
        'KRIGspline2','KRIGsphere2','KRIGcub2','POLYquad'})) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) <2*Data.dim+1+(Data.dim-1)*Data.dim/2) 
        fprintf('With Kriging and second order polynomial I need at least\n 2*Data.dim+1+(Data.dim-1)*Data.dim/2 starting points.\n')
        number_startpoints =2*Data.dim+1+(Data.dim-1)*Data.dim/2-size(starting_point,1);
    elseif(number_startpoints <2*Data.dim+1+(Data.dim-1)*Data.dim/2)  
        fprintf('With Kriging and second order polynomial I need at least\n 2*Data.dim+1+(Data.dim-1)*Data.dim/2 starting points.\n')
        number_startpoints =2*Data.dim+1+(Data.dim-1)*Data.dim/2;
    end
end

%check if there are enough starting points for using linear regression 
%model
if any(strcmp(surrogate_model,'POLYlin')) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) < 1+Data.dim)
        fprintf('For linear polynomial I need at least 1+Data.dim starting points.\n')
        number_startpoints =1+Data.dim-size(starting_point,1);
    elseif (number_startpoints <1+Data.dim)
        fprintf('For linear polynomial I need at least 1+Data.dim starting points.\n')
        number_startpoints =1+Data.dim;
    end
end

%check if there are enough starting points for using MARS model
if any(strcmp(surrogate_model,'MARS')) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) < 1+Data.dim)
        fprintf('For the MARS model I need at least 1+Data.dim starting points.\n')
        number_startpoints =1+Data.dim-size(starting_point,1);
    elseif (number_startpoints <1+Data.dim)
        fprintf('For the MARS model I need at least 1+Data.dim starting points.\n')
        number_startpoints =1+Data.dim;
    end
end

%check if there are enough starting points for using full cubic regression 
%model
if any(strcmp(surrogate_model,'POLYcub')) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) < 1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3))
        fprintf(['For full cubic polynomial I need at least\n 1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3) starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n']);
        number_startpoints =1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)-size(starting_point,1);
    elseif (number_startpoints <1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)) 
        fprintf(['For full cubic polynomial I need at least\n 1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3) starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n']);
        number_startpoints =1+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3);
    end
end

%check if there are enough starting points for using reduced cubic 
%polynomial regression model
if any(strcmp(surrogate_model,'POLYcubr')) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) < 1+3*Data.dim)
        fprintf(['For reduced cubic polynomial I need at least 1+3*Data.dim starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'])
        number_startpoints =1+3*Data.dim-size(starting_point,1);
    elseif (number_startpoints <1+3*Data.dim) 
        fprintf(['For reduced cubic polynomial I need at least 1+3*Data.dim starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n']);
        number_startpoints =1+3*Data.dim;
    end
end

%check if there are enough starting points for using reduced quadratic 
%polynomial regression model
if any(strcmp(surrogate_model,'POLYquadr')) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) < 1+2*Data.dim)
        fprintf(['For reduced quadratic polynomial I need at least 1+2*Data.dim starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'])
        number_startpoints =1+2*Data.dim-size(starting_point,1);
    elseif (number_startpoints <1+2*Data.dim) 
        fprintf(['For reduced quadratic polynomial I need at least 1+2*Data.dim starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'])
        number_startpoints =1+2*Data.dim;
    end
end

%check if there are enough starting points for using mixtures with full
%cubic polynomial regression model
if any(strcmp(surrogate_model,{'MIX_RcPc','MIX_KgPc'})) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) < 2+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)) %2+ because is mixture model
        fprintf(['For mixture with full cubic polynomial I need at least\n 2+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3) starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'])
        number_startpoints =2+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)-size(starting_point,1);
    elseif (Data.dim >=3)&& (number_startpoints <2+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3)) 
        fprintf(['For mixture with full cubic polynomial I need at least\n 2+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3) starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'])
        number_startpoints =2+3*Data.dim+ nchoosek(Data.dim,2)+ nchoosek(Data.dim,3);
    elseif Data.dim <3
        fprintf('For 2-dim problem there are no intersection terms of the form x1x2x3 possible.\n');
    end
end

%check if there are enough starting points for using mixtures with reduced
%quadratic polynomial regression model
if any(strcmp(surrogate_model,{'MIX_KgPqr'})) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) <2*Data.dim+2) 
        fprintf(['With Kriging and second order polynomial I need at least\n 2*Data.dim+1+(Data.dim-1)*Data.dim/2 starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'])
        number_startpoints =2*Data.dim+2-size(starting_point,1);
    elseif(number_startpoints <2*Data.dim+2)  
        fprintf(['With Kriging and second order polynomial I need at least\n 2*Data.dim+1+(Data.dim-1)*Data.dim/2 starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'])
        number_startpoints =2*Data.dim+2;
    end
end

%check if there are enough starting points for using mixtures with reduced
%cubic polynomial regression model
if any(strcmp(surrogate_model,{'MIX_KgPcr'})) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) < 2+3*Data.dim) %2+ because is mixture model
        fprintf(['For mixture with reduced cubic polynomial I need at least\n 2+3*Data.dim starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'])
        number_startpoints =2+3*Data.dim-size(starting_point,1);
    elseif (number_startpoints <2+3*Data.dim) 
        fprintf(['For mixture with reduced cubic polynomial I need at least\n 2+3*Data.dim starting points.\n',...
            'Otherwise the least squares problem is underdetermined.\n'])
        number_startpoints =2+3*Data.dim;
    end
end

%check if there are enough starting points for using mixtures with MARS or
%radial basis functions
if any(strcmp(surrogate_model,{'MIX_RcKg','MIX_RcM','MIX_KgM','MIX_RcKgM'})) 
    if exist('starting_point','var') && (number_startpoints+size(starting_point,1) < Data.dim+2) %2+ because is mixture model
        fprintf('For the chosen mixture I need at least 2+Data.dim starting points.\n')
        number_startpoints =2+Data.dim-size(starting_point,1);
    elseif (number_startpoints <2+Data.dim) 
        fprintf('For the chosen mixture I need at least 2+Data.dim starting points.\n')
        number_startpoints =2+Data.dim;
    end
end

%check if there are enough starting points for using SLHD design
if strcmp(initial_design,'SLHD') && exist('starting_point','var')
    if (number_startpoints+size(starting_point,1) <2*Data.dim)
        fprintf('With SLHD I need to use at least 2*Data.dim starting points.\n')
        number_startpoints = 2*Data.dim-size(starting_point,1);
    end
elseif strcmp(initial_design,'SLHD') 
    if (number_startpoints < 2*Data.dim)
        fprintf('With SLHD I need to use at least 2*Data.dim starting points.\n')
        number_startpoints = 2*Data.dim;
    end
    
end

%check if problem has integrality constraints. If so, make sure candidate
%point sampling is used
if exist('sampling_technique','var') && ~isempty(sampling_technique)
    if isempty(Data.continuous)|| (~isempty(Data.continuous) && ~isempty(Data.integer)) %pure or mixed-integer problem
        if ~strcmp(sampling_technique,'CAND')
            fprintf('For problems with integrality constraints you can only use the candidate\n point sampling.\n')
            sampling_technique='CAND';
        end
    end
end


%check if any input arguments were left empty and assign default values if
%needed
%also check if desired surrogate model is implemented
if isempty(surrogate_model) %surrogate model
    fprintf('Using default surrogate model.\n')
    surrogate_model='RBFcub';   %if no surrogate model defined, use cubic RBF
else %add own surrogate model to the list
     if ~any(strcmp(surrogate_model,{'RBFcub', 'RBFtps', 'RBFlin','KRIGexp0','KRIGexp1','KRIGexp2',...
             'KRIGgexp0', 'KRIGgexp1', 'KRIGgexp2','KRIGgauss0', 'KRIGgauss1', 'KRIGgauss2',...
             'KRIGlin0','KRIGlin1','KRIGlin2','KRIGspline0','KRIGspline1','KRIGspline2',...
             'KRIGsphere0','KRIGsphere1','KRIGsphere2','KRIGcub0','KRIGcub1','KRIGcub2',...
             'POLYlin', 'POLYquad', 'POLYquadr', 'POLYcub', 'POLYcubr', 'MARS', 'MIX_RcKg',...
             'MIX_RcM','MIX_RcPc', 'MIX_KgM', 'MIX_KgPc', 'MIX_KgPqr', 'MIX_KgPcr',...
             'MIX_RcKgM'}))
         error('The surrogate model you want to use is not contained in the toolbox. Check spelling.')
     end
end
if isempty(sampling_technique) %sampling technique
    fprintf('Using default sampling strategy.\n')
    sampling_technique='CAND';  %if no sampling method definded, use candidate points
else
    if ~any(strcmp(sampling_technique,{'CAND' 'SURFmin', 'EImax',...
            'SCOREmin', 'BUMPmin'}))
        error('The sampling technique you want to use is not contained in the toolbox. Check spelling.')
    end  
end
if isempty(initial_design) %initial experimental design
    fprintf('Using default initial experimental design strategy.\n')
    initial_design='LHS';      %if no initial experimental design strategy defined, use symm. Latin hypercube
else
    if ~any(strcmp(initial_design,{'SLHD', 'LHS', 'CORNER', 'SPACEFIL'}))
        error('The experimentala design strategy you want to use is not contained in the toolbox. Check spelling.')
    end
end
if isempty(maxeval) %maximum number of function evaluations
    fprintf('Using default number of maximal allowed function evaluations.\n')
    maxeval= 400;               %if maximum number of function evaluations not defined, set to 400
end
if isempty(number_startpoints) %starting point(s)
    fprintf('Using default number of points for initial experimental design.\n')
    number_startpoints=2*(Data.dim+1); %default is number of starting points not defined
end

%generate initial experimental design
Data.S=StartingDesign(initial_design,number_startpoints,Data);  %create initial experimental design
if exist ('starting_point', 'var') %user gives one or more starting points to add to the initital design
    Data.S=[starting_point; Data.S];
end

%regenerate initial experimental design if rank of sample site matrix too
%low
if isempty(Data.continuous) && ~isempty(Data.integer) %purely integer problem
    %round integer variables
    Data.S=round(Data.S);
    while rank([Data.S,ones(size(Data.S,1),1)]) <Data.dim+1 %regenerate starting design if rank of matrix too small
        Data.S=StartingDesign(initial_design,number_startpoints,Data);  %create initial experimental design
        if exist ('starting_point', 'var') %user gives one or more starting points to add to the initital design
            Data.S=[starting_point; Data.S];
        end
        Data.S=round(Data.S);
    end    
elseif ~isempty(Data.continuous) && isempty(Data.integer) %purely continuous problem 
    while rank([Data.S,ones(size(Data.S,1),1)]) <Data.dim+1 %regenerate starting design if rank of matrix too small
        Data.S=StartingDesign(initial_design,number_startpoints, Data);  %create initial experimental design
        if exist ('starting_point', 'var') %user gives one or more starting points to add to the initital design
            Data.S=[starting_point; Data.S];
        end
    end
elseif ~isempty(Data.continuous) && ~isempty(Data.integer) %mixed-integer problem    
    %round integer variables
    Data.S(:,Data.integer)=round(Data.S(:,Data.integer));
    while rank([Data.S,ones(size(Data.S,1),1)]) < Data.dim+1 %regenerate starting design if rank of matrix too small
        Data.S=StartingDesign(initial_design,number_startpoints, Data);  %create initial experimental design
        if exist ('starting_point', 'var') %user gives one or more starting points to add to the initital design
            Data.S=[starting_point; Data.S];
        end
        Data.S(:,Data.integer)=round(Data.S(:,Data.integer));
    end 
else
    error('Please enter variable vectors defining which variables are integers/continuous')
end

%for constrained mixed-integer problems, check if feasible starting point
%given
if ~isempty(Data.continuous) && ~isempty(Data.integer) %MI problem
    if isfield(Data,'constraint')  && (~exist('starting_point','var') || isempty(starting_point))
        fprintf('For a mixed-integer problem with black-box constraints\n you must provide at least one feasible starting point.\n')
        return
    elseif exist('starting_point','var') 
        if any (starting_point < Data.xlow) || any (starting_point >Data.xup)
            fprintf('The starting point is not within the box-constrained variable domain.\n')
            return
        end
        if isfield(Data,'constraint')     
            g=zeros(size(starting_point,1),length(Data.constraint));
            for kk=1:length(Data.constraint)
                g(:,kk)=feval(Data.constraint{kk},starting_point);
            end
            if ~any(sum(max(g,zeros(size(g))),2)<=0)
                fprintf('None of the starting points is feasible.\n')
            end
        end
    end
end

%do expensive function evaluations at points in initial experimental design
Data.Y=zeros(size(Data.S,1),1);
for ii = 1:size(Data.S,1)
    fevalst=tic; %record function evaluation time 
    Data.Y(ii,1)=feval(Data.objfunction,Data.S(ii,:));
    Data.fevaltime(ii,1)=toc(fevalst);
end 


%%----parallel version of doing function evalutaions, delete comments if
%%parallel version is wanted
% fvals=zeros(size(Data.S,1),1);
% object=Data.objfunction;
% Samples=Data.S;
% Tfeval=zeros(size(Data.S,1),1);
% %matlabpool(2)
% parfor ii = 1:size(Data.S,1)
%     fevalst=tic; %record function evaluation time 
%     fvals(ii,1)=feval(object,Samples(ii,:));
%     Tfeval(ii,1)=toc(fevalst);
% end  
% Data.Y=fvals;
% Data.fevaltime=Tfeval;
%--------end of parallel version-------------------

%optimization phase - distinguish between continuous/integer/mixed-integer
if isempty(Data.continuous) && ~isempty(Data.integer) %purely integer problem
    Data= OptimizationPhase_integer(Data,surrogate_model,maxeval);    
elseif ~isempty(Data.continuous) && isempty(Data.integer) %purely continuous problem 
    Data= OptimizationPhase_continuous(Data,surrogate_model,sampling_technique,maxeval);   
elseif ~isempty(Data.continuous) && ~isempty(Data.integer) %mixed-integer problem   
    Data= OptimizationPhase_mixedinteger(Data,surrogate_model,maxeval);    
end

Data.Problem=data_file;
Data.SurrogateModel=surrogate_model;
Data.SamplingTechnique=sampling_technique;
Data.InitialDesign=initial_design;
Data.NumberStartPoints=number_startpoints;
if exist('starting_point','var')
    Data.StartingPoint=starting_point;
end
Data.TotalTime=toc(timeStart);
save('Results','Data');

end %function