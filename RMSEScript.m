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
     InternalRMSEValues = 1:size(OutputArray(i).ys,2);
     InternalTrainingPointLength = 1:size(OutputArray(i).ys,2);
     for j = 1:size(OutputArray(i).ys,2)
        InternalYS = transpose(OutputArray(i).ys(:,j));
        InternalParam = OutputArray(i).parameters;
     
        TestPointLength = size(InternalParam(j).TestPoints,1);
        InternalYSEval = 1:TestPointLength;
        TestPointEval = 1:TestPointLength;
     
        for k = 1:TestPointLength
            TestPointEval(k) = feval(InternalParam(j).objectFunction, InternalParam(j).TestPoints(k,:));
        end
     
        InternalTrainingPointLength(j) = InternalParam.numStartPoints
        
        if(InternalParam(j).numStartPoints > 50)
            display("TESTING")
        end
        
        InternalRMSEValues(j) = getRMSE(InternalYS, TestPointEval);
     end
     RMSEValues = [RMSEValues InternalRMSEValues];
     TrainingPointLength = [TrainingPointLength InternalTrainingPointLength];
 end
 
 Figure1 = figure();
 subplot(2,1,1)
 plot(TrainingPointLength, real(RMSEValues), 'rx')
 xlabel("Training Point Length")
 ylabel("re(RMSE)")
 title("re(RMSE) vs. Training Point Length")
 subplot(2,1,2)
 plot(TrainingPointLength, imag(RMSEValues), 'bs')
 xlabel("Training Point Length")
 ylabel("im(RMSE)")
 title("im(RMSE) vs. Training Point Length")
