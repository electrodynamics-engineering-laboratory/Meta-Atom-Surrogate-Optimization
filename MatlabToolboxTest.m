clear;
%{
Inputs:
    Matlab Surrogate Model Module Test
    data_file: File name of function that establishes upper and lower
            limits, object function to be used by model, indices of integer
            variables, indices of continuous variables, and the problem dimension.
            See datainput_SOBModelProject.m for firmer example of syntax required.
    Num_Iterations: Number of iterations to be used by
                    SurrogateModelModule_v1.m
    SBOModels: Array of strings indicating the models to be tested with given
                inputs.
    Samp_Tech: String indicating the sampling technique to be used by model
    Init_Design: String indicating the initial design to be used by model
    Num_Start_Pnts: Integer indicating the number of starting points to be
                    used by the model
    Start_Point: An MxD matrix of points to be used by the model as initial
                    points where D is the dimension and M is an integer
                    greater than or equal to one.
Outputs:
  ResultOutput: An array of structures containing data output by each
                individual model

Joe Haun
Dr. Mohamed Salem
3-20-19

%}

%%
if(exist('ToolboxTestInputs.mat') ~= 2)
    data_file = ["datainput_SBOModel1" "datainput_SBOModel2"];
    Num_Iterations = 200;
    SBOModels = ["KRIGexp0" "KRIGexp1" "KRIGexp2" "KRIGgexp0" "KRIGgexp1" "KRIGgexp2" "KRIGgauss0" "KRIGgauss1" "KRIGgauss2" ...
        "KRIGlin0" "KRIGlin1" "KRIGlin2" "KRIGspline0" "KRIGspline1" "KRIGspline2" "KRIGsphere0" "KRIGsphere1" "KRIGsphere2" ...
        "KRIGcub0" "KRIGcub1" "KRIGcub2"];
    Samp_Tech = ["CAND", "SURFmin", "EImaxmin", "SCOREmin"];
    Init_Design = ["LHS", "SLHD", "CORNER", "SPACEFIL"];
    Num_Start_Pnts = 50;
    Start_Point = randn(10,4);
    save('ToolboxTestInputs.mat', 'data_file','Num_Iterations','SBOModels','Samp_Tech','Init_Design','Num_Start_Pnts','Start_Point');
else
    load('ToolboxTestInputs.mat');
end

%%
%Initialize ResultsOutput structure array and ErrorLog array
TestScalingFactor = 1;
StructureLengths = TestScalingFactor*length(Samp_Tech)*length(Init_Design)*length(SBOModels);
ResultOutput = struct([]);
ErrorLog = strings(2, StructureLengths);

%Expand ResultOutput array to prevent excessive slowdown when running tests
for i = 1:StructureLengths
    for j = 1:2
   ResultOutput(j,i).xlow = [];
   ResultOutput(j,i).xup = [];
   ResultOutput(j,i).objfunction = [];
   ResultOutput(j,i).integer = [];
   ResultOutput(j,i).continuous = [];
   ResultOutput(j,i).dim = [];
   ResultOutput(j,i).S = [];
   ResultOutput(j,i).Y = [];
   ResultOutput(j,i).fevaltime = [];
   ResultOutput(j,i).fbest = [];
   ResultOutput(j,i).xbest = [];
   ResultOutput(j,i).Ymed = [];
   ResultOutput(j,i).Problem = [];
   ResultOutput(j,i).SurrogateModel = [];
   ResultOutput(j,i).SamplingTechnique = [];
   ResultOutput(j,i).InitialDesign = [];
   ResultOutput(j,i).NumberStartPoints = [];
   ResultOutput(j,i).StartingPoint = [];
   ResultOutput(j,i).TotalTime = [];
    end
end
%%
AskForInput = 0; %If 1, check for user input before proceeding to next model
BreakOut = 0;
for ITERATIONS = 1:TestScalingFactor
    for k = 1:(length(Samp_Tech))
        for j = 1:(length(Init_Design))
            for i = 1:(length(SBOModels))
                for fileChoice = 1:length(data_file)
                    %display(SBOModels(i))
                    try %Enter try catch loop to prevent tests that fail from ending run of the program
                        SurrogateModelModule_v1(char(data_file(fileChoice)), Num_Iterations, char(SBOModels(i)), char(Samp_Tech(k)), char(Init_Design(j)), Num_Start_Pnts, Start_Point);
                        TempRes = load("Results.mat");
                        ResultOutput(fileChoice,ITERATIONS*i*j*k) = TempRes.Data;
                        ErrorLog(fileChoice,ITERATIONS*i*j*k) = strcat(string(Samp_Tech(k)), "-", string(Init_Design(j)), "-", string(SBOModels(i)), ": ", "OK");
                    catch ME
                        ErrorString = strcat(string(Samp_Tech(k)), "-", string(Init_Design(j)), "-", string(SBOModels(i)), ": ", string(ME.message));
                        ErrorLog(fileChoice,ITERATIONS*i*j*k) = ErrorString;
                    end
                end
            end
        end
    %{
    for k = 1:(length(Samp_Tech))
        for j = 1:(length(Init_Design))
            for i = 1:(length(SBOModels))
                %display(SBOModels(i))
                try %Enter try catch loop to prevent tests that fail from ending run of the program
                    SurrogateModelModule_v1(char(data_file(2)), Num_Iterations, char(SBOModels(i)), char(Samp_Tech(k)), char(Init_Design(j)), Num_Start_Pnts, Start_Point);
                    TempRes = load("Results.mat");
                    ResultOutput(2,ITERATIONS*i*j*k) = TempRes.Data;
                    ErrorLog(2,ITERATIONS*i*j*k) = strcat(string(Samp_Tech(k)), "-", string(Init_Design(j)), "-", string(SBOModels(i)), ": ", "OK");
                catch ME
                    ErrorString = strcat(string(Samp_Tech(k)), "-", string(Init_Design(j)), "-", string(SBOModels(i)), ": ", string(ME.message));
                    ErrorLog(2,ITERATIONS*i*j*k) = ErrorString;
                end
            end
        end
    end
    %}
    end
end
%%
OutFile = strcat(char(datetime),'_','ToolboxTestResults.mat');
save(OutFile, 'ResultOutput', 'ErrorLog' );
