function Ymixtest=model2combi2_3mod(Intervals,yPred_all)

%computes weights for models in 3-model mixtures and the predictions of the
%mixture model
%
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
%Intervals - matrix with believe, plausibilities, pignistic probabilities
%for each model in the mixture
%yPred_all - matrix with predictions of single surrogate models
%
%output:
%Ymixtest - mixture model prediction
%--------------------------------------------------------------------------

a1=find(Intervals(:,1)==1);
a2=find(Intervals(:,2)==2);
a3=find(Intervals(:,3)==3);

Inew=Intervals([a1,a2,a3],[1:3,end]);

Ymixtest=zeros(size(yPred_all,1),4); %(m x 11) matrix combinations in colns: 1/2 1/3 1/4 2/3 2/4 3/4 1/2/3 1/2/4 1/3/4 2/3/4 1/2/3/4
% contains predictions of mixed models at sampled sites
Mweight2=[];
for ii = 1:2
    for jj=ii+1:3
        weights=Inew([ii,jj],end)./sum(Inew([ii,jj],end));
        Mweight2=[Mweight2; ii jj weights(:)'];    %weights for combining 2 models   
    end
end
for ii = 1:size(Mweight2,1)
    Ymixtest(:,ii)=Mweight2(ii,3)*yPred_all(:,Mweight2(ii,1))+Mweight2(ii,4)*yPred_all(:,Mweight2(ii,2)); %predictions at sample points in leave one out
end

%mixture model prediction
Mweight3=Inew(:,end)'./sum(Inew(:,end));
Ymixtest(:,end)=Mweight3(1,1)*yPred_all(:,1)+Mweight3(1,2)*yPred_all(:,2)+Mweight3(1,3)*yPred_all(:,3);

end %function