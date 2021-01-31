function w=weights_in_combi(Intervals,newmodel,Data)
%determine weights for models in mixture
%--------------------------------------------------------------------------
%Copyright (c) 2012 by Juliane Mueller
%
% This file is part of the surrogate model module toolbox.
%
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%Tampere University of Technology, Finland
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%Intervals - matrix, contains mixture models with plausibilities, believes,
%and pignistic probabilities
%newmodel - vector, contains indices of the models in the mixture
%Data - structure, contains all model information
%
%output:
%w - vector with weights of each model contributing to combination
%
%--------------------------------------------------------------------------

I = NaN*ones(length(newmodel),Data.numberOfModels);
for ii=1:length(newmodel)
    I(ii,newmodel(ii))=newmodel(ii);
end

%determine weights of each model in the combined model based on
%pignistic probabilities
w=zeros(length(newmodel),1);
for ii = 1:length(newmodel)
    for jj = 3:size(Intervals,1)
        if  (sum(isnan(Intervals(jj,1:Data.numberOfModels)))==sum(isnan(I(ii,:)))) &&...
                (all(~isnan(Intervals(jj,1:Data.numberOfModels))-~isnan(I(ii,:))==0)) 
            w(ii,1)=Intervals(jj,end);
        end
    end
end
end %function