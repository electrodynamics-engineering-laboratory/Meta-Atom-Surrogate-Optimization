clear; close all;

%Get all file names to parse out .mat files containing data
FileStructure = dir("TestingOutputs");
FileStructureLength = length(FileStructure);
FileNames = string(zeros(1,FileStructureLength));
DesiredExtension = '.mat';
ExcludedFile = 'ToolboxTestInputs';
for i = 1:FileStructureLength
    [Path, Name, Extension] = fileparts(FileStructure(i).name);
    if(strcmp(Extension,DesiredExtension) && ~strcmp(Name,ExcludedFile))
        FileNames(i) = string(FileStructure(i).name); 
    end
end

FileNames = [FileNames(FileNames ~= "0")];
Figures = [];

for i = 1:length(FileNames)
    load(FileNames(i));
    
    OriginalErrorLogLength = length(ErrorLog); %Save initial length for calculations
    
    ErrorLog = rmmissing([ErrorLog(1,:)]); %Remove <missing> entries from log for parsing
    ErrorLog = ErrorLog(ErrorLog ~= ""); %Remove "" entries from log for parsing
    
    ZerosString = string(zeros(1,length(ErrorLog))); 
    ParsedLog.SamplingTechnique = ZerosString;
    ParsedLog.Models = ZerosString;
    ParsedLog.InitialDesign = ZerosString;
    ParsedLog.Status = ZerosString;
    clear ZerosString;

    for j = 1:length(ErrorLog)
        InternalTempString = char(ErrorLog(j));
        ColonLocations = find(InternalTempString == ':');
        ParsedLog.SamplingTechnique(j) = string(InternalTempString(1:(ColonLocations(1)-1)));
        ParsedLog.InitialDesign(j) = string(InternalTempString((ColonLocations(1)+1):(ColonLocations(2)-1)));
        ParsedLog.Models(j) = string(InternalTempString((ColonLocations(2)+1):(ColonLocations(3)-1)));
        ParsedLog.Status(j) = string(InternalTempString((ColonLocations(3)+1):end));
    end
    
    SuccessVal = 0;
    FailVal = 0;
    for k = 1:length(ErrorLog)
       SuccessVal = SuccessVal + strcmp(ParsedLog.Status(k), " OK");
    end
    FailVal = OriginalErrorLogLength - SuccessVal;
    
    Figures = [Figures figure('units','normalized','outerposition', [0 0 1 1])];
    subplot(2,2,1)
    bar(transpose([SuccessVal 0]));
    hold on
    bar([0 FailVal])
    hold off
    title("SBO Overall Succeses/Failures (Missing Data Included)");
    xticklabels(["Success", "Fail"]);
    legend(["Success", "Failure"]);
    pause(0.5)
    
    SBOModels = ["KRIGexp0" "KRIGexp1" "KRIGexp2" "KRIGgexp0" "KRIGgexp1" "KRIGgexp2" "KRIGgauss0" "KRIGgauss1" "KRIGgauss2" "KRIGlin0" "KRIGlin1" "KRIGlin2" "KRIGspline0" "KRIGspline1" "KRIGspline2" "KRIGsphere0" "KRIGsphere1" "KRIGsphere2" "KRIGcub0" "KRIGcub1" "KRIGcub2"];
    Samp_Tech = ["CAND", "SURFmin", "EImax", "SCOREmin"];
    Init_Design = ["LHS", "SLHD", "SPACEFIL"];

    ModelsSuccess = zeros(1,length(SBOModels));
    ModelsFailure = zeros(1,length(SBOModels));
    
    for j = 1:length(SBOModels)
        for k = 1:length(ErrorLog)
            if(strcmp(ParsedLog.Models(k), SBOModels(j)))
                ModelsSuccess(j) = ModelsSuccess(j) + strcmp(ParsedLog.Status(k), " OK");
                ModelsFailure(j) = ModelsFailure(j) + ~strcmp(ParsedLog.Status(k), " OK");
            end
        end
    end
    
    subplot(2,2,2)
    bar(transpose([ModelsSuccess; ModelsFailure]));
    title("SBO Models Success/Failure (Missing Data Excluded)");
    xticks(1:length(SBOModels));
    xticklabels(SBOModels);
    xtickangle(90);
    legend(["Success", "Failure"]);
    if(max(ModelsFailure) > max(ModelsSuccess))
        ylim([0 (max(ModelsFailure)+5)]);
    else
        ylim([0 (max(ModelsSuccess)+5)]);
    end
    pause(0.5);
    
    SampleSuccess = zeros(1,length(Samp_Tech));
    SampleFailure = zeros(1,length(Samp_Tech));
    
    for j = 1:length(Samp_Tech)
       for k = 1:length(ErrorLog)
            if(strcmp(ParsedLog.SamplingTechnique(k), Samp_Tech(j)))
                SampleSuccess(j) = SampleSuccess(j) + strcmp(ParsedLog.Status(k), " OK");
                SampleFailure(j) = SampleFailure(j) + ~strcmp(ParsedLog.Status(k), " OK");
            end
        end
    end
    
    subplot(2,2,3);
    bar(transpose([SampleSuccess; SampleFailure]));
    title("Sampling Technique Success/Failure (Missing Data Excluded)");
    xticks(1:length(Samp_Tech));
    xticklabels(Samp_Tech);
    xtickangle(90);
    legend(["Success", "Failure"]);
    if(max(SampleFailure) > max(SampleSuccess))
        ylim([0 (max(SampleFailure)+5)]);
    else
        ylim([0 (max(SampleSuccess)+5)]);
    end
    pause(0.5);
    
    InitialDesignSuccess = zeros(1,length(Init_Design));
    InitialDesignFailure = zeros(1,length(Init_Design));
    
    for j = 1:length(Init_Design)
        for k = 1:length(ErrorLog)
            if(strcmp(ParsedLog.InitialDesign(k), Init_Design(j)))
                InitialDesignSuccess(j) = InitialDesignSuccess(j) + strcmp(ParsedLog.Status(k), " OK");
                InitialDesignFailure(j) = InitialDesignFailure(j) + ~strcmp(ParsedLog.Status(k), " OK");
            end
        end
    end
    
    subplot(2,2,4)
    bar(transpose([InitialDesignSuccess; InitialDesignFailure]));
    title("Initial Design Success/Failure  (Missing Data Excluded)");
    xticks(1:length(Init_Design));
    xticklabels(Init_Design);
    xtickangle(90);
    legend(["Success", "Failure"]);
    if(max(InitialDesignFailure) > max(InitialDesignSuccess))
        ylim([0 (max(InitialDesignFailure)+5)]);
    else
        ylim([0 (max(InitialDesignSuccess)+5)]);
    end
    pause(0.5);

end


