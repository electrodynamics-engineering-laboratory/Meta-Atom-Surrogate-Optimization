function Output = runModelTest(inputStartPoints)
%{
    This function is for running tests on the EELMetamodel/get_ys_krig
    functions. 

    
%}

% Setup variables for running get_ys_krig function
if nargin == 0
    inputStartPoints = 1; %If no number is given, assume one test is to be run
end

for numStartPoints = (1):(inputStartPoints)  
    
    testsToRun = 1; %Only run the specific parameters once
    func = @(x)real((x(:,3)+i*x(:,4))./(x(:,1)+i*x(:,2) + x(:,3)+i*x(:,4))); %The function to be tested against
    dimension = 4; %The function's number of variables
    %numStartPoints = 40; %The number of start points which are used to create training points
    xLow = [0.2,0.0,0,0]; 
    xHigh = [10,10,10,10];
    thetaLow = 0.1;
    thetaHigh = 0.9;

    %Create arrays for all possible initial designs, sampling strategies, and regression polynomials
    initDesign = ["lhsdesign"];
    sampStrat = ["Cubic", "Exponential", "Gaussian", "Linear", "Spherical", "Spline"];
    regrPoly = [0, 1, 2];

    %Preallocate a blank structure for data input later
    Output = struct('ys', [], 'dmodel', [], 'error', [], 'parameters', []);

    %Preallocate dummy structures in case of error, this separates outputs
    %for later parsing. 
    tempDModel = struct('regr', NaN, 'corr', NaN, 'theta', NaN, 'beta', NaN, 'gamma', NaN, 'sigma2', NaN, 'S', NaN, 'Ssc', NaN, 'Ysc', NaN, 'C', NaN, 'Ft', NaN, 'G', NaN);
    tempParameters = struct('TrainingPoints', NaN, 'FunctionValues', NaN, 'TestPoints', NaN*zeros(4), 'dimension', NaN, 'numStartPoints', NaN, 'initialDesign', NaN, 'xLim', NaN, 'objectFunction', func, 'thetaBound', NaN, 'correlation', NaN, 'regressionPolynomial', NaN);
    tempYs = NaN*zeros(dimension,1);           
    
    for Iterations = 1:testsToRun
        for IDChoice = 1:length(initDesign)
           for SSChoice = 1:length(sampStrat)
              for RPChoice = 1:length(regrPoly)
                  try
                    tempStruct = struct('initDesign', initDesign(IDChoice), 'sampStrat', sampStrat(SSChoice), 'regrPoly', regrPoly(RPChoice), 'message', 'NONE');
                    [ys, dmodel, parameters] = EELMetamodel(dimension, numStartPoints, initDesign(IDChoice), xLow, xHigh, func, thetaLow, thetaHigh, sampStrat(SSChoice), regrPoly(RPChoice));
                    Output.ys = [Output.ys ys];
                    Output.dmodel = [Output.dmodel dmodel];
                    Output.parameters = [Output.parameters parameters];
                    Output.error = [Output.error tempStruct];
                  catch MatlabError
                    tempStruct = struct('initDesign', initDesign(IDChoice), 'sampStrat', sampStrat(SSChoice), 'regrPoly', regrPoly(RPChoice), 'message', string(MatlabError.message));
                    Output.ys = [Output.ys tempYs];
                    Output.parameters = [Output.parameters tempParameters];
                    Output.dmodel = [Output.dmodel tempDModel];
                    Output.error = [Output.error tempStruct];
                  end
              end %RPChoice Loop
           end %SSChoice Loop
        end %IDChoice Loop
    end %testsToRun Loop

% Save the output structure into a file
    OutputLocation = "../TestingOutputs/";
    TestName = strcat("DACE_MAT_",string(numStartPoints),"TP_");
    DateString = char(datetime);
    DateString(DateString == ' ') = '_';
    DateString(DateString == ':') = '-';
    OutFile = strcat(OutputLocation,TestName,DateString,'_','OutputResults.mat');
    save(OutFile, 'Output');

end %inputStartPoints Loop

end


