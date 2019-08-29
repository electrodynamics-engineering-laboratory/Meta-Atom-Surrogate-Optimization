function [AccuracyOutput] = evaluateAccuracy(EELMetamodelOutput)
%{
    This function takes in the output of the runModelTest function and
    creates graphs of the ys values and the object function evaluated at
    the test points

    Input: EELMetamodelOutput - an array of structures containing the
                                fields: ys, dmodel, error, and parameters
                    
                    ys - an m-by-n matrix of values where each column
                    corresponds to a dimension of the problem and each row
                    corresponds to a specific run of the get_ys_krig
                    function
                    dmodel - an array of structures that is the output of 
                    the dacefit function which contains various parameters 
                    used to create a metamodel
                    error - an array that contains any error messages
                    output by the runModelTest function
                    parameters - an array of structures that contains the
                    values input to the get_ys_krig function within the
                    EELMetamodel function
    
    Output: AccuracyOutput - an array of structures containing the fields:
                    TestPointEval, GivenYSValues, NormalizedYS, and Figures
                    
                    TestPointEval - a m-by-n matrix of values created by
                    evaluating the given object function at each testing
                    point given
                    GivenYSValues - an m-by-n matrix of values that are
                    given to the function
                    NormalizedYS - an m-by-n matrix of values that are each
                    element of the GivenYSValues divided by each element of
                    the TestPointEval
                    
%}

Input = EELMetamodelOutput;
AccuracyOutput = struct('TestPointEval', [], 'GivenYSValues', [], 'NormalizedYS', [], 'Figures', []);

for k = 1:size(Input.ys,2)

InputYS = transpose(Input.ys(:,k));
InputParameters = Input.parameters(:,k);

TestPointLength = size(InputParameters.TestPoints,1);
TestPointEval = 1:TestPointLength;

for i = 1:TestPointLength
    TestPointEval(i) = feval(InputParameters.objectFunction, InputParameters.TestPoints(i,:));
end

NormalizedYSValues = InputYS./TestPointEval;

AccuracyOutput.TestPointEval = [AccuracyOutput.TestPointEval; TestPointEval];
AccuracyOutput.GivenYSValues = [AccuracyOutput.GivenYSValues; InputYS];
AccuracyOutput.NormalizedYS = [AccuracyOutput.NormalizedYS; NormalizedYSValues];

end

for k = 1:size(Input.ys,1)
   TempFigure = figure();
   plot([1:length(AccuracyOutput.GivenYSValues(:,k))], AccuracyOutput.GivenYSValues(:,k), 'ro')
   hold on
   plot([1:length(AccuracyOutput.NormalizedYS(:,k))], AccuracyOutput.TestPointEval(:,k), 'bx')
   hold off
   legend('Given YS Values', 'Test Point Evaluations', 'Normalized YS Values', 'Location', 'southoutside')
   title(strcat("Dimension Parameter (",string(k),")"));
   xlabel("Evaluation Number");
   ylabel(strcat("Value"));
   AccuracyOutput.Figures = [AccuracyOutput.Figures TempFigure];
   pause(0.5)
   
end



end