function [RootMeanSquareError] = calculateRMSE(ActualValues, TheoreticalValues)
%{
    [RootMeanSquareError] = calculateRMSE(ActualValues, TheoreticalValues)
    
    Function takes in actual values output and theoretical values output
    then determines the Root Mean Square Error.

    Inputs:
        -ActualValues: An MxN matrix of the values observed/collected
        -TheoreticalValues: An MxN matrix of the values calculated
    Outputs:
        -RootMeanSquareError: 
%}
ActualLength = length(ActualValues);
TheoreticalLength = length(TheoreticalValues);

if(ActualLenth ~= TheoreticalLength)
    display("Error: Input Lengths Do Not Match")
    RootMeanSquareError = -1;
    return  
else
    RootMeanSquareError = sqrt((sum(ActualValues - TheoreticalValues).^2)/ActualLength);
end

end