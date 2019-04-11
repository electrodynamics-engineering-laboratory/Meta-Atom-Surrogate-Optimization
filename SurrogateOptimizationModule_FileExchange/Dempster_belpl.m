function IntervalsD=Dempster_belpl(BPA_dempster,Data)
%calculates believes, plausibilities and pignistic probabilities of focal
%elements based on Dempster's rule of combination
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
%BPA_dempster - matrix with basic probability assignments
%Data - structure with all model information
%Output:
%IntervalsD - matrix with belive, plausibility and pignistic probability
%values
%newmodelD - vector with models in mixture
%--------------------------------------------------------------------------

b=(1:Data.numberOfModels);

%% believes 
belD=zeros(size(BPA_dempster,1),1);
for ii = 1:size(BPA_dempster,1)
    cur1=b(~isnan(BPA_dempster(ii,1:Data.numberOfModels)));
    for jj = 1:size(BPA_dempster,1)
        cur2=b(~isnan(BPA_dempster(jj,1:Data.numberOfModels)));
        if all(ismember(cur2,cur1))
            belD(ii)=belD(ii)+BPA_dempster(jj,end);
        end
    end
end

%% plausibilities
plD=zeros(size(BPA_dempster,1),1);
for ii = 1:size(BPA_dempster,1)
    cur1=b(~isnan(BPA_dempster(ii,1:Data.numberOfModels)));
    for jj = 1:size(BPA_dempster,1)
        cur2=b(~isnan(BPA_dempster(jj,1:Data.numberOfModels)));
        if ~isempty(intersect(cur2,cur1))
            plD(ii)=plD(ii)+BPA_dempster(jj,end);
        end
    end
end
    
%% Pignistic probabilities (see FLOREA JOUSSELME BOSSE GRENIER: robust combination rules for evidence theory)
betP_D=zeros(size(BPA_dempster,1),1);
for ii = 1:size(BPA_dempster,1)
    cur1=b(~isnan(BPA_dempster(ii,1:Data.numberOfModels)));
    for jj = 1:size(BPA_dempster,1)
        cur2=b(~isnan(BPA_dempster(jj,1:Data.numberOfModels)));
        if ~isempty(intersect(cur2,cur1))
            betP_D(ii)=betP_D(ii)+BPA_dempster(jj,end)*length(intersect(cur2,cur1))/length(cur2);
        end
    end
end

%% believe -plausibility Intervals
IntervalsD=[BPA_dempster(:,1:Data.numberOfModels) belD plD betP_D];
IntervalsD=sortrows(IntervalsD,-(Data.numberOfModels+2));

end %function