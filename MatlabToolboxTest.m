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
    SOBModels: Array of strings indicating the models to be tested with given
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
Dr. Salem
3-20-19

%}

if(exist('ToolboxTestInputs.mat') ~= 2)
    data_file = 'datainput_SBOModelProject';
    Num_Iterations = 200;
    SOBModels = ["KRIGexp0" "KRIGexp1" "KRIGexp2" "KRIGgexp0" "KRIGgexp1" "KRIGgexp2" "KRIGgauss0" "KRIGgauss1" "KRIGgauss2" ...
        "KRIGlin0" "KRIGlin1" "KRIGlin2" "KRIGspline0" "KRIGspline1" "KRIGspline2" "KRIGsphere0" "KRIGsphere1" "KRIGsphere2" ...
        "KRIGcub0" "KRIGcub1" "KRIGcub2"];
    Samp_Tech = ["CAND", "SURFmin", "EImaxmin", "SCOREmin"];
    Init_Design = ["LHS", "SLHD", "CORNER", "SPACEFIL"];
    Num_Start_Pnts = 50;
    Start_Point = randn(10,4);
    save('ToolboxTestInputs.mat', 'data_file','Num_Iterations','SOBModels','Samp_Tech','Init_Design','Num_Start_Pnts','Start_Point');
else
    load('ToolboxTestInputs.mat');
end

ResultOutput = [];
ErrorLog = [];

AskForInput = 0; %If 1, check for user input before proceeding to next model
BreakOut = 0;
for k = 1:(length(Samp_Tech))
    for j = 1:(length(Init_Design))
        for i = 1:(length(SOBModels))
            display(SOBModels(i))
            try %Enter try catch loop to prevent tests that fail from ending run of the program
                SurrogateModelModule_v1(data_file, Num_Iterations, char(SOBModels(i)), char(Samp_Tech(k)), char(Init_Design(j)), Num_Start_Pnts, Start_Point);
                TempRes = load("Results.mat");
                ResultOutput = [ResultOutput TempRes.Data];
                ErrorString = strcat(string(Samp_Tech(k)), "-", string(Init_Design(j)), "-", string(SOBModels(i)), ": OK");
                ErrorLog = [ErrorLog ErrorString];
            catch ME
                ErrorString = strcat(string(Samp_Tech(k)), "-", string(Init_Design(j)), "-", string(SOBModels(i)), ": ", string(ME.message));
                ErrorLog = [ErrorLog ErrorString];
            end
            
            if(AskForInput == 1 && string(input('Continue? ', 's')) == "N") %If not continuing set variable to break out of loops
                BreakOut = 1;
                break;
            end
        end
        if(BreakOut == 1)
            break;
        end
    end
    if(BreakOut == 1)
        break;
    end
end

save('ToolboxTestResults.mat', 'ResultOutput', 'ErrorLog' );
