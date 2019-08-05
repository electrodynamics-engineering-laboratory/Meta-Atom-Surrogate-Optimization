function [Output] = evaluateAccuracy(EELMetamodelOutput)
%{

%}

Input = EELMetamodelOutput;
Output = struct('TestPointEval', [], 'GivenYSValues', [], 'NormalizedYS', [], 'Figures', []);

for k = 1:size(Input.ys,2)

InputYS = transpose(Input.ys(:,k));
InputParameters = Input.parameters(:,k);

TestPointLength = size(InputParameters.TestPoints,1);
TestPointEval = 1:TestPointLength;

for i = 1:TestPointLength
    TestPointEval(i) = feval(InputParameters.objectFunction, InputParameters.TestPoints(i,:));
end

NormalizedYSValues = InputYS./TestPointEval;

Output.TestPointEval = [Output.TestPointEval; TestPointEval];
Output.GivenYSValues = [Output.GivenYSValues; InputYS];
Output.NormalizedYS = [Output.NormalizedYS; NormalizedYSValues];

end

for k = 1:size(Input.ys,1)
   TempFigure = figure();
   plot([1:length(Output.GivenYSValues(:,k))], Output.GivenYSValues(:,k), 'ro')
   hold on
   plot([1:length(Output.NormalizedYS(:,k))], Output.TestPointEval(:,k), 'bx')
   hold off
   legend('Given YS Values', 'Test Point Evaluations', 'Location', 'southoutside')
   title(strcat("X",string(k)));
   xlabel("Evaluation Number");
   ylabel(strcat("Value"));
   Output.Figures = [Output.Figures TempFigure];
   
end

end