function Figures = ParseResultsOutput(Data, GivenFileName)
%{
    Function takes in an array of structures from the
    SurrogateModelToolboxTest function and generates appropriate figures to
    display the data gathered

    Inputs: Data - an array of structures output by the
    SurrogateModelToolboxTest function
    Outputs: Figures - an array of figure handles corresponding to all
    generated figures

    Author: Joe Haun
    6-14-19

%}
    %Establish Function Constants
    fbestString = 'fbest';
    xbestString = 'xbest';
    totalTimeString = 'TotalTime';
    if nargin < 2
       GivenFileName = "UNKNOWN"; 
    end
    Figures = [];
    DataFileNames = fieldnames(Data);
    TempArray = [];
    for i = 1:length(DataFileNames)
       TempArray = [TempArray string(DataFileNames{i})];
    end
    DataFileNames = TempArray;
    
    for i = 1:length(DataFileNames)
       eval(strcat("CurDataFileFields = fieldnames(Data.",DataFileNames(i),");"));
       ParsedData = struct(fbestString, [], xbestString, [], totalTimeString, []);
       for j = 1:length(CurDataFileFields)
           %Find desired fields
           if strcmp(CurDataFileFields{j}, fbestString)
               eval(strcat("ParsedData.", fbestString," = [Data.", DataFileNames(i), ".",fbestString,"];"));
           elseif strcmp(CurDataFileFields{j}, xbestString)
               eval(strcat("TempArray = 1:length([Data.",DataFileNames(i),"]);"));
               if strcmp(DataFileNames(i), "datainput_Schoen_10_4_3")
                  DataFileNames(i);
               end
               for k = TempArray
                  eval(strcat("ParsedData.",xbestString, " = [ParsedData.", xbestString, " transpose([Data.", DataFileNames(i), "(k).", CurDataFileFields{j}, "])];"))
                  k
               end
           elseif strcmp(CurDataFileFields{j}, totalTimeString)
               eval(strcat("ParsedData.", totalTimeString, " = [Data.", DataFileNames(i), ".", totalTimeString, "];"));
           end
       end
       %Generate new figure for subplots
       if( ~isempty(ParsedData.TotalTime) || ~isempty(ParsedData.fbest) || ~isempty(ParsedData.xbest))
            Figures = [Figures figure('Name', GivenFileName ,'units','normalized','outerposition', [0 0 1 1])];
       end
       
       %Generate subplots of data
       if( ~isempty(ParsedData.TotalTime))
            subplot(2,2,1);
            plot([ParsedData.TotalTime], 'o')
            xlabel("Run Number")
            ylabel("Time (sec)")
            title(strcat(DataFileNames(i), " Total Time"));
            pause(0.25);
       end
       
       if( ~isempty(ParsedData.fbest))
            subplot(2,2,2);
            plot(real(ParsedData.fbest), imag(ParsedData.fbest), 'o')
            xlabel("fbest (Real)")
            ylabel("fbest (Imag)")
            title(strcat(DataFileNames(i), " fBest"));
            pause(0.25);
       end
       
       if( ~isempty(ParsedData.xbest))
            subplot(2,2,[3, 4]);
            plot(1:length(ParsedData.xbest(1,:)), ParsedData.xbest(1,:), 'o', 1:length(ParsedData.xbest(2,:)), ParsedData.xbest(2,:),  'x', 1:length(ParsedData.xbest(3,:)), ParsedData.xbest(3,:), '+', 1:length(ParsedData.xbest(4,:)), ParsedData.xbest(4,:), 's') 
            xlabel("")
            ylabel("Time (sec)")
            title(strcat(DataFileNames(i), " xBest"));
            legend(["R1", "X1", "R2", "X2"]);
            pause(0.25);
       end
       
    end
       
       
end %End of Function

