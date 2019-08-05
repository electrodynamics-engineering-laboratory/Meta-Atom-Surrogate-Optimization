function [Figures] = evaluateAccuracy(EELMetamodelOutput)
%{

%}

Input = EELMetamodelOutput;

for k = 1:size(Input.ys,2)

InputYS = Input.ys(:,k);
InputParameters = Input.parameters(:,k);

TestPointLength = size(InputParameters.TestPoints,1);
TestPointEval = 1:TestPointLength;

for i = 1:TestPointLength
    TestPointEval(i) = feval(InputParameters.objectFunction, InputParameters.TestPoints(i,:));
end

TestPointEval

end


end