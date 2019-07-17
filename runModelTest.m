function Output = runModelTest(testsToRun)
%{

%}

%%
% Setup variables for running get_ys_krig function
if nargin == 0
    testsToRun = 1; %If no number is given, assume one test is to be run
end
func = @(x)((x(:,3)+i*x(:,4))./(x(:,1)+i*x(:,2) + x(:,3)+i*x(:,4))); %The function to be tested against
dimension = 4; %The function's number of variables
numStrPnt = 10; %The number of start points for the function? (Might be wrong)
xLow = [0.1,0.1,0.1,0.1]; 
xHigh = [10,10,10,10];
thetaLow = 0.1;
thetaHigh = 0.9;

%Create arrays for all possible initial designs, sampling strategies, and regression polynomials
initDesign = ["lhsdesign"];
sampStrat = ["Cubic", "Exponential", "Gaussian", "Linear", "Spherical", "Spline"];
regrPoly = [0, 1, 2];

%Preallocate a blank structure for data input later
Output = struct('ys', [], 'dmodel', [], 'error', []);

%%
for Iterations = 1:testsToRun
    for IDChoice = 1:length(initDesign)
       for SSChoice = 1:length(sampStrat)
          for RPChoice = 1:length(regrPoly)
              try
                [ys, dmodel] = EELMetamodel(dimension, numStrPnt, initDesign(IDChoice), xLow, xHigh, func, thetaLow, thetaHigh, sampStrat(SSChoice), regrPoly(RPChoice));
                Output.ys = [Output.ys ys];
                Output.dmodel = [Output.dmodel dmodel];
              catch MatlabError
                tempStruct = struct('initDesign', initDesign(IDChoice), 'sampStrat', sampStrat(SSChoice), 'regrPoly', regrPoly(RPChoice), 'message', string(MatlabError.message));
                Output.error = [Output.error tempStruct];
              end
          end
       end
    end
end

%%
% Save the output structure into a file
OutputLocation = "TestingOutputs/";
TestName = "DACE";
DateString = char(datetime);
DateString(DateString == ' ') = '_';
DateString(DateString == ':') = '-';
OutFile = strcat(OutputLocation,TestName,DateString,'_','OutputResults.mat');
save(OutFile, 'Output');

end


