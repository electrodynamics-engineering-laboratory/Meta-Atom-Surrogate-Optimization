function Figures = DisplayResults(FileLocation)
%{
    Function takes in an directory and then creates graphics of both the 
    Results and Errors reported by the SurrogateModelToolboxTest function
    
    Inputs:     OutputFileLocation - A string of the directory path
    Outputs:    Figures - An array of figure handles for any graphs created
%}
DesiredExtension = ".mat";
Figures = [];
[filepath, name, extension] = fileparts(string(FileLocaton));
if strcmp(extension,DesiredExtension)
    load(FileLocation); %Load current file as it meets desired criteria
    Figures = [Figures ParseResultsOutput(ResultOutput,name) ParseErrorLog(ErrorLog, name)]; %Append generated figures to figures array
end

end