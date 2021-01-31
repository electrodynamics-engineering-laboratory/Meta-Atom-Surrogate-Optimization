function modelweights = DempsterFor3models(Data, Surrogate)
%Dempster-Shafer-Theory for combining 3 models
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
%Data - structure with all problem information
%Surrogate - string with name of mixture surrogate model
%Output:
%modelweights - vector with weights for each model
%--------------------------------------------------------------------------

m=size(Data.S,1); %numer of already sampled points
yPred_all=zeros(m,3); %initialize matrix for predicted function values

%% determine whether k-fold or leave-one-out cross-validation
if m>50  %k-fold cross validation
    multiplier=ceil(m/50)-1;
    kfold=10*multiplier;
    forcekfold=true;
    mix=randperm(m);
else %leave-one-out cross-validation
    forcekfold=false;
end 

Data.numberOfModels=3; %number of models in mixture
if strcmp(Surrogate,'MIX_RcKgM')
    if ~forcekfold  %leave-one-out cross-validation
        for jj=1:m
            %RBF model prediction
            [lambda, gamma]=RBF([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cubic'); 
            yPred_all(jj,1)=RBF_eval(Data.S(jj,:),[Data.S(1:jj-1,:);...
                Data.S(jj+1:end,:)],lambda,gamma,'cubic');
            %Kriging model prediction
            dmodel = dacefit([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)], 'regpoly1', ...
                'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
            yPred_all(jj,2)=predictor(Data.S(jj,:), dmodel);  
            %MARS model prediction
            mmodel = aresbuild([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)]); 
            yPred_all(jj,3)=arespredict(mmodel,Data.S(jj,:)); 
        end      
    else
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);
   
        for jj=1:ceil(m/kfold)
            %validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            % training set
            if jj ==1 %first group
                trainset_S=Ss(jj*kfold+1:end,:);
                trainset_Y=Ys(jj*kfold+1:end,:);
            elseif 2<=jj && jj<= floor(m/kfold)  && mod(m,kfold) >0 ||...
                    2<=jj && mod(m,kfold)==0 && jj*kfold<m
                    %group 2 till last full group %mod(m,kfold)~=0
                trainset_S=[Ss(1:(jj-1)*kfold,:);Ss(jj*kfold+1:end,:)];
                trainset_Y=[Ys(1:(jj-1)*kfold,:);Ys(jj*kfold+1:end,:)];
            else %left-overs
                trainset_S=Ss(1:(jj-1)*kfold,:);
                trainset_Y=Ys(1:(jj-1)*kfold,:);
            end
            if 1<=jj && jj<= floor(m/kfold)  && mod(m,kfold) >0 ||...
                    2<=jj && mod(m,kfold)==0 && jj*kfold<m
                %RBF model predictions
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                yPred_all((jj-1)*kfold+1:jj*kfold,1)=RBF_eval(validation_S,...
                    trainset_S,lambda,gamma,'cubic');
                %Kriging model predictions
                dmodel = dacefit(trainset_S, trainset_Y, 'regpoly1', 'corrgauss',...
                    ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=predictor(validation_S,dmodel);
                %MARS model predictions
                mmodel = aresbuild(trainset_S,trainset_Y); 
                yPred_all((jj-1)*kfold+1:jj*kfold,3)=arespredict(mmodel,validation_S);
            else %left-overs
                id_a=(jj-1)*kfold+1;
                if mod(m,kfold)>0
                    id_e=(jj-1)*kfold+mod(m,kfold);
                else
                    id_e=jj*kfold;
                end
                %RBF model predictions
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                yPred_all(id_a:id_e,1)=RBF_eval(validation_S,trainset_S,...
                    lambda,gamma,'cubic');
                %Kriging model predictions
                dmodel=dacefit(trainset_S, trainset_Y, 'regpoly1', 'corrgauss',...
                	ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                yPred_all(id_a:id_e,2)=predictor(validation_S,dmodel);
                %MARS model predictions
                mmodel = aresbuild(trainset_S,trainset_Y); 
                yPred_all(id_a:id_e,3)=arespredict(mmodel,validation_S);
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end  
    end
%add further 3-model combinations here     
%elseif strcmp(Surrogate,'YourModel')
end
        
%% characteristics of individual models 
Corcoef=cc_calc(Data.Ymed,yPred_all); %correlation coefficient
RMSE=RMSE_calc(Data.Ymed,yPred_all); %root mean squared error
MAE=MAE_cal(Data.Ymed,yPred_all); %maximum absolute errors
MAD=MAD_cal(Data.Ymed,yPred_all); %median absolute deviation
    
%correlation coefficients    
Corcoef_t=Corcoef;
if any(Corcoef_t<0) %scale correlation coefficients to interval 0,1
    Corcoef_t=(Corcoef_t-min(Corcoef_t))./max(Corcoef_t-min(Corcoef_t));
end
if all(isnan(Corcoef))
    Corcoef_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    Corcoef_t(Corcoef_t==0)=0.000001;
    Corcoef_t=Corcoef_t(:)./(sum(Corcoef_t(:)));
end
%Root mean squared errors
if all(RMSE==0)
    RMSE_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    RMSE_t=(1./RMSE(:))./sum(1./RMSE(:));
end
% Maximum absolute errors
if all(MAE==0)
    MAE_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    MAE_t=(1./MAE(:))./sum(1./MAE(:));
end
%median absolute deviation
if all(MAD==0)
    MAD_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    MAD_t=(1./MAD(:))./(sum(1./MAD));
end

%combine model information
comb=[1 NaN NaN; NaN 2 NaN; NaN NaN 3]; %individual models
OVpcr=[comb Corcoef_t(:), RMSE_t(:), MAE_t(:), MAD_t(:)]; %evidences of individual models
    
intround=true;
%DST to combine evidences
for yy =2:Data.numberOfModels
    if intround
        compD=OVpcr(:,1:Data.numberOfModels+1);
        intround=false;
    else
        compD=Intersections;
    end
    OVnew=[comb,OVpcr(:,Data.numberOfModels+yy)];%Transform data into evidences
    Intersections=dempster_rule(compD,OVnew,Data);%combination rule
end
%calculate believes, plausibilities and pignistic probabilities
Intervals=Dempster_belpl(Intersections,Data);
%predict function values with all possible 2-model combinations and with 3-model combination    
Ymixtest=model2combi2_3mod(Intervals,yPred_all); 
%characteristics of mixture
Corcoef=[Corcoef;cc_calc(Data.Ymed,Ymixtest)];
if any(Corcoef<0)
    Corcoef=(Corcoef-min(Corcoef))./max(Corcoef-min(Corcoef));
end
Corcoef(Corcoef==0)=0.000001;
RMSE=[RMSE;RMSE_calc(Data.Ymed,Ymixtest)]; %root mean squared errors
MAE=[MAE;MAE_cal(Data.Ymed,Ymixtest)]; %maximum absolute errors
MAD=[MAD;MAD_cal(Data.Ymed,Ymixtest)]; % median absolute deviation


%% Combine pure models to mixed model     
%all available model combinations
comb=[1 NaN NaN; NaN 2 NaN; NaN NaN 3; 1 2 NaN; 1 NaN 3; NaN 2 3; 1 2 3]; 

%matrix containing model combination and corresponding evidences (high numbers mean high evidence)
OVpcr=[comb Corcoef(:)./(sum(Corcoef(:))) (1./RMSE(:))./sum(1./RMSE(:))...
    (1./MAE(:))./sum(1./MAE(:)) (1./MAD(:))./(sum(1./MAD))];
clear Corcoef; clear RMSE; clear MAE; clear MAD
clear Corcoef_t; clear RMSE_t; clear MAE_t; clear MAD_t;

intround=true;
%DST
for yy =2:Data.numberOfModels
    if intround
        compD=OVpcr(:,1:Data.numberOfModels+1);
        intround=false;
    else
        compD=Intersections;
    end
    OVnew=[comb,OVpcr(:,Data.numberOfModels+yy)]; %Transform data into evidences
    Intersections=dempster_rule(compD,OVnew,Data); %combination rule
end

%calculate believes, plausibilities and pignistic probabilities
Intervals=Dempster_belpl(Intersections,Data);
newmodel=Intervals(1,1:Data.numberOfModels);
newmodel(isnan(newmodel))=[];
%determine the weights of the models in the combination
w=weights_in_combi(Intervals,newmodel,Data);
modelweights=zeros(1,Data.numberOfModels);
modelweights(newmodel)=w/sum(w);
end%function