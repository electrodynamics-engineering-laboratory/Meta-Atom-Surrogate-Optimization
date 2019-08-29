clear;
TargetPath = "TestingOutputs";
TargetDirectory = dir(TargetPath);
ValidFiles = string(1:length(TargetDirectory));

for i = 1:length(TargetDirectory)
    if(~strcmp(TargetDirectory(i).name, '.') & ~strcmp(TargetDirectory(i).name,'..') & ~strcmp(TargetDirectory(i).name, 'Old_Runs'))
        ValidFiles(i) = strcat(TargetPath, "/", string(TargetDirectory(i).name));
    else
       ValidFiles(i) = "000"; 
    end
end

OutputArray = struct('ys', [], 'dmodel', [], 'error', [], 'parameters', []);

 for i = 1:length(ValidFiles)
     try
        load(ValidFiles(i)); 
        OutputArray = [OutputArray Output];
     catch ME
         display(ME.message)
     end
 end
 
 %Calculate RMSE Values and Graph
 RMSEValues = [];
 TrainingPointLength = [];
 
 %Get TrainingPoint lengths
 for i = 2:length(OutputArray)
     InternalRMSEValues = 1:length(OutputArray);
     InternalTrainingPointLength = 1:length(OutputArray);
     for j = 1:size(OutputArray(i).ys,2)
        InternalYS = transpose(OutputArray(i).ys(:,j));
        InternalParam = OutputArray(i).parameters;
     
        TestPointLength = size(InternalParam(j).TestPoints,1);
        InternalYSEval = 1:TestPointLength;
        TestPointEval = 1:TestPointLength;
     
        for k = 1:TestPointLength
            TestPointEval(k) = feval(InternalParam(j).objectFunction, InternalParam(j).TestPoints(k,:));
        end
     
        InternalTrainingPointLength(j) = InternalParam.numStartPoints;
        InternalRMSEValues(j) = getRMSE(InternalYS, TestPointEval);
     end
     RMSEValues = [RMSEValues transpose(InternalRMSEValues)];
     TrainingPointLength = [TrainingPointLength transpose(InternalTrainingPointLength)];
 end