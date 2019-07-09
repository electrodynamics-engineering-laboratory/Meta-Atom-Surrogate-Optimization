function Output = runModelTest()
%Run EELMetamodel
%Attempts to run script with following command
dimension = 4;
numStrPnt = 10;
xLow = 0.1;
xHigh = 1;
func = @(x)((x(:,3)+i*x(:,4))./(x(:,1)+i*x(:,2) + x(:,3)+i*x(:,4)));
thetaLow = 0.1;
thetaHigh = 0.9;

initDesign = ["lhsdesign"];
sampStrat = ["Cubic", "Exponential", "Gaussian", "Linear", "Spherical", "Spline"];
regrPoly = [0, 1, 2];

Output = struct('ys', [], 'dmodel', [], 'error', []);
Output.dmodel = struct('regr', [], 'corr', [], 'theta', [], 'beta', [], 'gamma', [], 'sigma2', [], 'S', [], 'Ssc', [], 'Ysc', [], 'C', [], 'Ft', [], 'G', []);
Output.error = struct('initDesign', [], 'sampStrat', [], 'regrPoly', [], 'message', []);
for IDChoice = 1:length(initDesign)
   for SSChoice = 1:length(sampStrat)
      for RPChoice = 1:length(regrPoly)
          try
            [ys, dmodel] = EELMetamodel(dimension, numStrPnt, initDesign(IDChoice), xLow, xHigh, func, thetaLow, thetaHigh, sampStrat(SSChoice), regrPoly(RPChoice));
            Output.ys = [Output.ys ys];
            Output.dmodel = [Output.dmodel dmodel];
          catch MatlabError
              tempStruct = struct('initDesign', initDesign(IDChoice), 'sampStrat', sampStrat(SSChoice), 'regrPoly', regrPoly(RPChoice), 'message', string(MatlabError.message));
            Output.error = [Output.error tempStruct];
          end
      end
   end
end



end


