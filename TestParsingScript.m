function Figures = ParseErrorLogs(OutputFileLocation)
%{
    Function takes in a directory and then creates graphics based on the
    error logs located within those .mat files to display relevant
    information in regards to the runs done.

    Inputs:     OutputFileLocation - A string of the directory path
    Outputs:    Figures - An array of figure handles for any graphs created
%}
if nargin < 1
    OutputFileLocation = "TestingOutputs";
end

%Get all file names to parse out .mat files containing data
FileStructure = dir(OutputFileLocation);
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
idealFontSize = 18;

for i = length(FileNames):length(FileNames)
    load(FileNames(i));

    ErrorLogFieldNames = fieldnames(ErrorLog);
    Temp = string(zeros(1,length(ErrorLogFieldNames)));
    for j = 1:length(ErrorLogFieldNames)
        Temp(j) = string(ErrorLogFieldNames{j});
    end

    ErrorLogLength = length(eval(strcat("ErrorLog.", ErrorLogFieldNames(1)))); %Save initial length for calculations

    ZerosString = string(zeros(1,ErrorLogLength)); 
    ParsedLog.SamplingTechnique = ZerosString;
    ParsedLog.Models = ZerosString;
    ParsedLog.InitialDesign = ZerosString;
    ParsedLog.Status = ZerosString;
    ParsedLog.FileName = ZerosString;
    clear ZerosString;

    ErrorLogFieldNames = Temp;
    clear Temp;

    for j = 1:length(ErrorLogFieldNames)
        TempArray = eval(strcat("[ErrorLog.",ErrorLogFieldNames(j),"]"));
        for k = 1:length(TempArray)
           InternalTempString = char(TempArray(k));
           ColonLocations = find(InternalTempString == ':');
           ParsedLog.SamplingTechnique(k) = string(InternalTempString(1:(ColonLocations(1)-1)));
           ParsedLog.InitialDesign(k) = string(InternalTempString((ColonLocations(1)+1):ColonLocations(2)-1));
           ParsedLog.Models(k) = string(InternalTempString((ColonLocations(2)+1):(ColonLocations(3)-1)));
           ParsedLog.Status(k) = string(InternalTempString((ColonLocations(3)+1):end));
           ParsedLog.FileName(k) = ErrorLogFieldNames(j);
        end
    end

    SuccessVal = 0;
    FailVal = 0;
    for k = 1:length(ParsedLog.Status)
       SuccessVal = SuccessVal + (strcmp(ParsedLog.Status(k), "OK") || strcmp(ParsedLog.Status(k), " OK")) ;
       FailVal = FailVal + ~(strcmp(ParsedLog.Status(k),"OK") || strcmp(ParsedLog.Status(k), " OK"));
    end

    Figures = [Figures figure('units','normalized','outerposition', [0 0 1 1])];
    title(FileNames(i));
    subplot(2,2,1)
    bar(transpose([SuccessVal 0]));
    hold on
    bar([0 FailVal])
    hold off
    title("SBO Overall Succeses/Failures");
    xticklabels(["Success", "Fail"]);
    AxesLegend = legend(["Success", "Failure"]);
    AxesLegend.FontSize = 12;
    curAxes = gca;
    curAxes.FontSize = idealFontSize;
    pause(0.25)

    SBOModels = ["KRIGexp0" "KRIGexp1" "KRIGexp2" "KRIGgexp0" "KRIGgexp1" "KRIGgexp2" "KRIGgauss0" "KRIGgauss1" "KRIGgauss2" "KRIGlin0" "KRIGlin1" "KRIGlin2" "KRIGspline0" "KRIGspline1" "KRIGspline2" "KRIGsphere0" "KRIGsphere1" "KRIGsphere2" "KRIGcub0" "KRIGcub1" "KRIGcub2"];
    Samp_Tech = ["CAND", "SURFmin", "EImax", "SCOREmin"];
    Init_Design = ["LHS", "SLHD", "SPACEFIL"];

    ModelsSuccess = zeros(1,length(SBOModels));
    ModelsFailure = zeros(1,length(SBOModels));

    for j = 1:length(SBOModels)
        for k = 1:ErrorLogLength
            if(strcmp(ParsedLog.Models(k), SBOModels(j)))
                ModelsSuccess(j) = ModelsSuccess(j) + (strcmp(ParsedLog.Status(k), "OK") || strcmp(ParsedLog.Status(k), " OK"));
                ModelsFailure(j) = ModelsFailure(j) + ~(strcmp(ParsedLog.Status(k), "OK") || strcmp(ParsedLog.Status(k), " OK"));
            end
        end
    end

    subplot(2,2,2)
    bar(transpose([ModelsSuccess; ModelsFailure]));
    title("SBO Models Success/Failure");
    xticks(1:length(SBOModels));
    xticklabels(SBOModels);
    xtickangle(90);
    AxesLegend = legend(["Success", "Failure"]);
    AxesLegend.FontSize = 12;
    if(max(ModelsFailure) > max(ModelsSuccess))
        ylim([0 (max(ModelsFailure)+max(ModelsFailure)/2)]);
    else
        ylim([0 (max(ModelsSuccess)+max(ModelsSuccess)/2)]);
    end
    curAxes = gca;
    curAxes.FontSize = idealFontSize
    pause(0.25);

    SampleSuccess = zeros(1,length(Samp_Tech));
    SampleFailure = zeros(1,length(Samp_Tech));

    for j = 1:length(Samp_Tech)
       for k = 1:ErrorLogLength
            if(strcmp(ParsedLog.SamplingTechnique(k), Samp_Tech(j)))
                SampleSuccess(j) = SampleSuccess(j) + (strcmp(ParsedLog.Status(k), "OK") || strcmp(ParsedLog.Status(k), " OK"));
                SampleFailure(j) = SampleFailure(j) + ~(strcmp(ParsedLog.Status(k), "OK") || strcmp(ParsedLog.Status(k), " OK"));
            end
        end
    end

    subplot(2,2,3);
    bar(transpose([SampleSuccess; SampleFailure]));
    title("Sampling Technique Success/Failure");
    xticks(1:length(Samp_Tech));
    xticklabels(Samp_Tech);
    AxesLegend = legend(["Success", "Failure"]);
    AxesLegend.FontSize = 12;
    if(max(SampleFailure) > max(SampleSuccess))
        ylim([0 (max(SampleFailure)+max(SampleFailure)/2)]);
    else
        ylim([0 (max(SampleSuccess)+max(SampleSuccess)/2)]);
    end
    curAxes = gca;
    curAxes.FontSize = idealFontSize
    pause(0.25);

    InitialDesignSuccess = zeros(1,length(Init_Design));
    InitialDesignFailure = zeros(1,length(Init_Design));

    for j = 1:length(Init_Design)
        for k = 1:ErrorLogLength
            if(strcmp(ParsedLog.InitialDesign(k), Init_Design(j)))
                InitialDesignSuccess(j) = InitialDesignSuccess(j) + (strcmp(ParsedLog.Status(k), "OK") || strcmp(ParsedLog.Status(k), " OK"));
                InitialDesignFailure(j) = InitialDesignFailure(j) + ~(strcmp(ParsedLog.Status(k), "OK") || strcmp(ParsedLog.Status(k), " OK"));
            end
        end
    end

    subplot(2,2,4)
    bar(transpose([InitialDesignSuccess; InitialDesignFailure]));
    title("Initial Design Success/Failure");
    xticks(1:length(Init_Design));
    xticklabels(Init_Design);
    AxesLegend = legend(["Success", "Failure"]);
    AxesLegend.FontSize = 12;
    if(max(InitialDesignFailure) > max(InitialDesignSuccess))
        ylim([0 (max(InitialDesignFailure)+max(InitialDesignFailure)/2)]);
    else
        ylim([0 (max(InitialDesignSuccess)+max(InitialDesignSuccess)/2)]);
    end
    curAxes = gca;
    curAxes.FontSize = idealFontSize
    pause(0.25);

end

end %End Function