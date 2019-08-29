function [RMSE] = getRMSE(ActualValues, TheoreticalValues)
%{
    [RootMeanSquareError] = getRMSE(ActualValues, TheoreticalValues)
    
    Function takes in actual values output and theoretical values output
    then determines the Root Mean Square Error and plots resulting values.

    Inputs:
        -ActualValues: An MxN matrix of the values observed/collected
        -TheoreticalValues: An MxN matrix of the values calculated
    Outputs:
        -RootMeanSquareError: 
%}
ActualLength = length(ActualValues);
TheoreticalLength = length(TheoreticalValues);

if(ActualLength ~= TheoreticalLength)
    display("Input Lengths Do Not Match")
    return  
end

RMSE = sqrt((sum(ActualValues - TheoreticalValues).^2)/ActualLength);

end