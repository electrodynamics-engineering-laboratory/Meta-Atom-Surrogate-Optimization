function Figures = DisplayResults(OutputFileLocation)
%{
    Function takes in an directory and then creates graphics of both the 
    Results and Errors reported by the SurrogateModelToolboxTest function
    
    Inputs:     OutputFileLocation - A string of the directory path
    Outputs:    Figures - An array of figure handles for any graphs created
%}
if nargin < 1
    OutputFileLocation = "TestingOutputs/";
end

Figures = [];

DesiredExtension = ".mat";
DirectoryList = dir(OutputFileLocation);
for i = 1:length(DirectoryList)
    [filepath, name, extension] = fileparts(string(DirectoryList(i).name));
    if strcmp(extension,DesiredExtension)
        load(strcat(OutputFileLocation,DirectoryList(i).name)); %Load current file as it meets desired criteria
        Figures = [Figures ParseResultsOutput(ResultOutput,DirectoryList(i).name) ParseErrorLog(ErrorLog, DirectoryList(i).name)]; %Append generated figures to figures array

    end
end


end