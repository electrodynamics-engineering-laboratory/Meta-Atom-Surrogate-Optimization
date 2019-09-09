function [ys, dmodel, Parameters] = EELMetamodel(dimension, numStartPoints, initialDesign, xLow, xHigh, objectFunction, thetaLowerBound, thetaUpperBound, correlation, regressionPolynomial)
%{
    Function creates the TrainingPoints, FunctionValues, and TestPoints for
    a Kriging Metamodel. The function accepts a structure with
    initial_design, dimension, ... as fields. 

    Input:
            dimension -            an integer value representing the total
                                   dimension of the problem (e.g. 1, 10, 32)
            numStartPoints -       an integer representing the total number
                                   of start points for the desired metamodel
            initialDesign -        a string of characters matching 'lhsdesign', 'SLHD',
                                   'bestlh', or 'cornerpoints'. Not case sensitive.
            xLow -                 an array of floating point values that
                                   indicate the lower bound for each 
                                   function parameter.
            xHigh -                an array of floating poitn values that
                                   indicate the upper bound for each
                                   function parameter.
            objectFunction  -      the function with which to compare the
                                   model and creating training points
            thetaLowerBound -      the lower limit of the theta value
            thetaUpperBound -      the upper limit of the theta value
            correlation -          Can be 'Cubic', 'Exponential', 'Gaussian', 'Linear', 'Spherical', or 'Spline'. Not case sensitive.  
            regressionPolynomial - An integer of either 0, 1, or 2
                                   corresponding to a 0th, 1st, or 2nd order.
                

    Output: 
                
%}

% Determine the initial design and create a matrix of training points depending on that choice
initialDesignChoices = ["lhsdesign",  "bestlh", "SLHD", "cornerpoints"];
if strcmpi( initialDesign,  initialDesignChoices(1))
    % Use the lhsdesign function built into matlab to provide trainig points
    TrainingPoints = lhsdesign(numStartPoints, dimension, 'criterion', 'maximin', 'iterations', 20);
    TrainingPoints = repmat(xHigh-xLow,numStartPoints,1).*TrainingPoints + repmat(xLow,numStartPoints,1);
elseif strcmpi( initialDesign,  initialDesignChoices(2))
    % Use the optimized latin hypercube to provide training points
    TrainingPoints = bestlh(numStartPoints,dimension,20,10);
    TrainingPoints = repmat(xHigh-xLow,size(TrainingPoints,1),1).*TrainingPoints + repmat(xLow,size(TrainingPoints,1),1);
elseif strcmpi( initialDesign,  initialDesignChoices(3))
    % Use the symmetric latin hypercube design to provide training points
    TrainingPoints = SLHD(dimension, numStartPoints);
    TrainingPoints = repmat(xHigh-xLow,size(TrainingPoints,1),1).*TrainingPoints + repmat(xLow,size(TrainingPoints,1),1);
elseif strcmpi( initialDesign,  initialDesignChoices(4))
    % Use the corner points design to provide training points
    TrainingPoints = cornerpoints(Conditions, numStartPoints);
else
    error("A valid initial design was not given. Aborting operation.")
    return;
end

% Set the correlation value to a number for use in the get_ys_krig function
correlationChoices = ["Cubic", "Exponential", "Gaussian", "Linear", "Spherical", "Spline"];
if strcmpi(correlation,correlationChoices(1))
    correlation = 1;
elseif strcmpi(correlation,correlationChoices(2))
    correlation = 2;
elseif strcmpi(correlation,correlationChoices(3))
    correlation = 3;
elseif strcmpi(correlation,correlationChoices(4))
    correlation = 4;
elseif strcmpi(correlation,correlationChoices(5))
    correlation = 5;    
elseif strcmpi(correlation,correlationChoices(6))
    correlation = 6;
else
    error("A valid correlation value was not given. Aborting operation.")
    return;
end

% Add one to the regressionPolynomial value to correspond with get_ys_krig function's desired inputs
if regressionPolynomial < 3 && regressionPolynomial >= 0
    regressionPolynomial = regressionPolynomial + 1;
else
    error("A valid regression Polynomial was not given. Aborting operation.")
    return;
end

% Given the amount of training points, evaluate the given function with each set of training points
FunctionValues = zeros(size(TrainingPoints,1),1);

for i = 1:size(TrainingPoints,1)
    FunctionValues(i,1) = feval(objectFunction, TrainingPoints(i,:));
end

% Create a set of random test points, this should be changed later to something better reflecting user input
%TestPoints = rand(dimension);
TestPoints = [1 2 3 4; 2 3 4 1; 3 4 1 2; 4 1 2 3];

[ys, dmodel] = get_ys_krig(TrainingPoints, FunctionValues, TestPoints, dimension, thetaLowerBound, thetaUpperBound, correlation, regressionPolynomial);

%Place all given parameters into the Parameters structure for output
Parameters = struct();
Parameters.TrainingPoints = TrainingPoints;
Parameters.FunctionValues = FunctionValues; 
Parameters.TestPoints = TestPoints;
Parameters.dimension = dimension;
Parameters.numStartPoints = numStartPoints;
Parameters.initialDesign = initialDesign;
Parameters.xLim = [xLow; xHigh];
Parameters.objectFunction = objectFunction;
Parameters.thetaBound = [thetaLowerBound thetaUpperBound];
Parameters.correlation = correlation;
Parameters.regressionPolynomial = regressionPolynomial;

end