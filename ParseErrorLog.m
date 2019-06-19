function Figures = ParseErrorLog(Data, GivenFileName)
%{
    Function takes in ErrorLog from run and parses out details

    Returns a structure with fields for each specific parameter
%}

if nargin < 2
   GivenFileName = "UNKNOWN"; 
end

%Get all file names to parse out .mat files containing data
idealFontSize = 18;
Figures = [];
DataFileNames = fieldnames(Data);
TempArray = [];
for i = 1:length(DataFileNames)
    TempArray = [TempArray string(DataFileNames{i})];
end
DataFileNames = TempArray;
TempArray = [];
for i = 1:length(DataFileNames)

    eval(strcat("TempArray = [Data.",DataFileNames{i},"];"))
    ZerosString = string(zeros(1,length(TempArray)));
    ParsedLog = struct('Models', ZerosString, 'SamplingTechnique', ZerosString, 'InitialDesign', ZerosString, 'Status', ZerosString, 'FileName', ZerosString);
    for j = 1:length(TempArray)

        InternalTempString = char(TempArray(j));
        ColonLocations = find(InternalTempString == ':');
        ParsedLog.SamplingTechnique(j) = string(InternalTempString(1:(ColonLocations(1)-1)));
        ParsedLog.InitialDesign(j) = string(InternalTempString((ColonLocations(1)+1):ColonLocations(2)-1));
        ParsedLog.Models(j) = string(InternalTempString((ColonLocations(2)+1):(ColonLocations(3)-1)));
        ParsedLog.Status(j) = string(InternalTempString((ColonLocations(3)+1):end));
        ParsedLog.FileName(j) = DataFileNames(i);

    end

    %Gather total successes and failures for current data file
    SuccessVal = 0;
    FailVal = 0;
    for j = 1:length(ParsedLog.Status)
       SuccessVal = SuccessVal + (strcmp(ParsedLog.Status(j), "OK") || strcmp(ParsedLog.Status(j), " OK")) ;
       FailVal = FailVal + ~(strcmp(ParsedLog.Status(j),"OK") || strcmp(ParsedLog.Status(j), " OK"));
    end

    %Gather successes and failures for each individual model
    SBOModels = ["KRIGexp0" "KRIGexp1" "KRIGexp2" "KRIGgexp0" "KRIGgexp1" "KRIGgexp2" "KRIGgauss0" "KRIGgauss1" "KRIGgauss2" "KRIGlin0" "KRIGlin1" "KRIGlin2" "KRIGspline0" "KRIGspline1" "KRIGspline2" "KRIGsphere0" "KRIGsphere1" "KRIGsphere2" "KRIGcub0" "KRIGcub1" "KRIGcub2"];
    Samp_Tech = ["CAND", "SURFmin", "EImax", "SCOREmin"];
    Init_Design = ["LHS", "SLHD", "SPACEFIL"];

    ModelsSuccess = zeros(1,length(SBOModels));
    ModelsFailure = zeros(1,length(SBOModels));

    for k = 1:length(SBOModels)
        for j = 1:length([ParsedLog.Status])
            if(strcmp(ParsedLog.Models(j), SBOModels(k)))
                ModelsSuccess(k) = ModelsSuccess(k) + (strcmp(ParsedLog.Status(j), "OK") || strcmp(ParsedLog.Status(j), " OK"));
                ModelsFailure(k) = ModelsFailure(k) + ~(strcmp(ParsedLog.Status(j), "OK") || strcmp(ParsedLog.Status(j), " OK"));
            end
        end
    end

    SampleSuccess = zeros(1,length(Samp_Tech));
    SampleFailure = zeros(1,length(Samp_Tech));

    for k = 1:length(Samp_Tech)
       for j = 1:length([ParsedLog.Status])
            if(strcmp(ParsedLog.SamplingTechnique(j), Samp_Tech(k)))
                SampleSuccess(k) = SampleSuccess(k) + (strcmp(ParsedLog.Status(j), "OK") || strcmp(ParsedLog.Status(j), " OK"));
                SampleFailure(k) = SampleFailure(k) + ~(strcmp(ParsedLog.Status(j), "OK") || strcmp(ParsedLog.Status(j), " OK"));
            end
        end
    end

    InitialDesignSuccess = zeros(1,length(Init_Design));
    InitialDesignFailure = zeros(1,length(Init_Design));

    for k = 1:length(Init_Design)
        for j = 1:length([ParsedLog.Status])
            if(strcmp(ParsedLog.InitialDesign(j), Init_Design(k)))
                InitialDesignSuccess(k) = InitialDesignSuccess(k) + (strcmp(ParsedLog.Status(j), "OK") || strcmp(ParsedLog.Status(j), " OK"));
                InitialDesignFailure(k) = InitialDesignFailure(k) + ~(strcmp(ParsedLog.Status(j), "OK") || strcmp(ParsedLog.Status(j), " OK"));
            end
        end
    end

    Figures = [Figures figure('Name', GivenFileName,'units','normalized','outerposition', [0 0 1 1])];
    title(DataFileNames(i));
    subplot(2,2,1)
    bar(transpose([SuccessVal 0]));
    hold on
    bar([0 FailVal])
    hold off
    title(strcat(DataFileNames{i},"SBO Overall Succeses/Failures"));
    xticklabels(["Success", "Fail"]);
    AxesLegend = legend(["Success", "Failure"]);
    AxesLegend.FontSize = 12;
    curAxes = gca;
    curAxes.FontSize = idealFontSize;
    pause(0.25)

    subplot(2,2,2)
    bar(transpose([ModelsSuccess; ModelsFailure]));
    title(strcat(DataFileNames{i},"SBO Models Success/Failure"));
    xticks(1:length(SBOModels));
    xticklabels(SBOModels);
    xtickangle(90);
    AxesLegend = legend(["Success", "Failure"]);
    AxesLegend.FontSize = 12;
    if(max(ModelsFailure) > max(ModelsSuccess))
        ylim([0 (max(ModelsFailure)+max(ModelsFailure)/2)]);
    else
        ylim([0 (1 + max(ModelsSuccess)+max(ModelsSuccess)/2)]);
    end
    curAxes = gca;
    curAxes.FontSize = idealFontSize;
    pause(0.25);

    subplot(2,2,3);
    bar(transpose([SampleSuccess; SampleFailure]));
    title(strcat(DataFileNames{i},"Sampling Technique Success/Failure"));
    xticks(1:length(Samp_Tech));
    xticklabels(Samp_Tech);
    AxesLegend = legend(["Success", "Failure"]);
    AxesLegend.FontSize = 12;
    if(max(SampleFailure) > max(SampleSuccess))
        ylim([0 (max(SampleFailure)+max(SampleFailure)/2)]);
    else
        ylim([0 (max(SampleSuccess)+max(SampleSuccess)/2 +1)]);
    end
    curAxes = gca;
    curAxes.FontSize = idealFontSize;
    pause(0.25);

    subplot(2,2,4)
    bar(transpose([InitialDesignSuccess; InitialDesignFailure]));
    title(strcat(DataFileNames{i},"Initial Design Success/Failure"));
    xticks(1:length(Init_Design));
    xticklabels(Init_Design);
    AxesLegend = legend(["Success", "Failure"]);
    AxesLegend.FontSize = 12;
    if(max(InitialDesignFailure) > max(InitialDesignSuccess))
        ylim([0 (max(InitialDesignFailure)+max(InitialDesignFailure)/2)]);
    else
        ylim([0 (1 + max(InitialDesignSuccess)+max(InitialDesignSuccess)/2)]);
    end
    curAxes = gca;
    curAxes.FontSize = idealFontSize;
    pause(0.25);

    end

end

