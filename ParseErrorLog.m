function ParsedErrors = ParseErrorLog(ErrorLog)
%{
    Function takes in ErrorLog from run and parses out details

    Returns a structure with fields for each specific parameter
%}

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

end

