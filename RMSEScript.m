clear;
TargetPath = "../TestingOutputs/";
TargetDirectory = dir(TargetPath);
ValidFiles = string(1:(length(TargetDirectory)-3));

minStringIndex = 10;
maxStringIndex = 13;
for i = 1:length(TargetDirectory)
    for j = 0:(maxStringIndex-minStringIndex)
        if(~isfolder(strcat(TargetPath,TargetDirectory(i).name)))    
            tempName = TargetDirectory(i).name(minStringIndex:(maxStringIndex-j));
            Position = find(ValidFiles == string(tempName));
            if(~isempty(Position))
               ValidFiles(Position) = strcat(TargetPath,TargetDirectory(i).name);
               break; 
            end
        end
    end
end

RMSEValues = zeros(length(ValidFiles), 18);
TrainingPointLengths = zeros(length(ValidFiles),18);
Correlation = zeros(length(ValidFiles),18);
RegressionPolynomial = zeros(length(ValidFiles),18);

%Get TrainingPoint length
 for i = 1:length(ValidFiles)
     load(ValidFiles(i))
     InternalRMSEValues = zeros(1, size(Output.ys,2));
     InternalTrainingPointLength = zeros(1,size(Output.ys,2));
     InternalCorrelation = zeros(1,size(Output.ys,2));
     InternalRegrPolynomial = zeros(1,size(Output.ys,2));
     
     for j = 1:size(Output.ys,2)
        InternalYS = transpose(Output.ys(:,j));
        InternalParam = Output.parameters;
     
        TestPointLength = size(InternalParam(j).TestPoints,1);
        InternalYSEval = 1:TestPointLength;
        TestPointEval = 1:TestPointLength;
     
        for k = 1:TestPointLength
            InternalParam(j).TestPoints(k,:);
            TestPointEval(k) = feval(InternalParam(j).objectFunction, InternalParam(j).TestPoints(k,:));
        end
        InternalCorrelation(j) = InternalParam(j).correlation;
        InternalRegrPolynomial(j) = InternalParam(j).regressionPolynomial;
        InternalTrainingPointLength(j) = InternalParam(j).numStartPoints;
        InternalRMSEValues(j) = getRMSE(InternalYS, TestPointEval);
     end
     RMSEValues(i,:) = InternalRMSEValues;
     TrainingPointLengths(i,:) = InternalTrainingPointLength;
     RegressionPolynomial(i,:) = InternalRegrPolynomial;
     Correlation(i,:) = InternalCorrelation;
 end
 ParamPairs = ["Cubic-Regr0", "Exponential-Regr0", "Gaussian-Regr0", "Linear-Regr0", "Spherical-Regr0", "Spline-Regr0" ...
       "Cubic-Regr1", "Exponential-Regr1", "Gaussian-Regr1", "Linear-Regr1", "Spherical-Regr1", "Spline-Regr1" ...
       "Cubic-Regr2", "Exponential-Regr2", "Gaussian-Regr2", "Linear-Regr2", "Spherical-Regr2", "Spline-Regr2"];
 Figures = [];
 for i = 1:size(RMSEValues,2)
     Figure = [Figures figure()];
     %subplot(2,1,1)
     plot(TrainingPointLengths(:,i), RMSEValues(:,i), 'rx')
     ylim([0, 0.25])
     xlim([0, (length(TargetDirectory)-3)])
     xlabel("Training Point Length")
     ylabel("re(RMSE)")
     title(strcat("re(RMSE) vs. Training Point Length: ", ParamPairs(i)))
 end