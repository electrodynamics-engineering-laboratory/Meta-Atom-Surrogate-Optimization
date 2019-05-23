clear;

%Get all file names to parse out .mat files containing data
FileStructure = dir;
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
    ErrorLog = [ErrorLog(ErrorLog ~= "")];
    ErrorLogCartesian = rmmissing([ErrorLog(1,:)]);
    ErrorLogPolar = rmmissing([ErrorLog(2,:)]);
    ErrorLog = [ErrorLogCartesian; ErrorLogPolar];
    keyboard()
    Figures = [Figures figure()];
    
    clear ResultOutput ErrorLog;
end



