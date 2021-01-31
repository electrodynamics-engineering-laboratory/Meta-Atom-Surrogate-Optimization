function modelweights = DempsterFor2models(Data,Surrogate)

%uses Dempster-Shafer Theory to determin the weight for model combinations
%with 2 models
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
%input: 
%Data - structure  with all model information so far
%Surrogate - desired mixture Surrogate model
%Output:
%modelweights - vector with weight of each model
%--------------------------------------------------------------------------

m=size(Data.S,1); %number of points already sampled
yPred_all=zeros(m,2);  %initialize array for objective function value predictions

%% determine whether k-fold or leave-one-out cross-validation
if m>50  %k-fold cross validation
    multiplier=ceil(m/50)-1;
    kfold=10*multiplier;
    forcekfold=true;
    mix=randperm(m);
else %leave-one-out cross-validation
    forcekfold=false;
end    

Data.numberOfModels=2;

%% cross-validation for mixtures
%mixture of cubic RBF and Kriging with Gaussian correlation
if strcmp(Surrogate,'MIX_RcKg')
    if ~forcekfold  %no k-fold, leave-1-out crossvalidation
        for jj=1:m %leave out data point jj
            %RBF model
            [lambda, gamma]=RBF([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cubic'); % RBF
            yPred_all(jj,1)=RBF_eval(Data.S(jj,:),[Data.S(1:jj-1,:);...
                Data.S(jj+1:end,:)],lambda,gamma,'cubic'); %repredict jj-th objective function value
            %Kriging model
            dmodel = dacefit([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)], 'regpoly1', ...
                'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  %Kriging
            if ~isstruct(dmodel)
                disp('Building the Kriging model failed')
                return
            else
                yPred_all(jj,2)=predictor(Data.S(jj,:), dmodel);   %repredict jj-th objective function value
            end
        end
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %randomly arranged sample sites
        Ys=Data.Ymed(mix);
        for jj=1:ceil(m/kfold)
            %determine points in validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            %determine points in training set
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
                % RBF model for prediction
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                yPred_all((jj-1)*kfold+1:jj*kfold,1)=RBF_eval(validation_S,...
                    trainset_S,lambda,gamma,'cubic');
                % Kriging model for prediction
                dmodel = dacefit(trainset_S, trainset_Y, 'regpoly1', 'corrgauss',...
                    ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                 yPred_all((jj-1)*kfold+1:jj*kfold,2)=predictor(validation_S,dmodel);
            else %for remaining points (left-overs)
                id_a=(jj-1)*kfold+1;
                if mod(m,kfold)>0
                    id_e=(jj-1)*kfold+mod(m,kfold);
                else
                    id_e=jj*kfold;
                end
                % RBF model for prediction
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                yPred_all(id_a:id_e,1)=RBF_eval(validation_S,trainset_S,...
                    lambda,gamma,'cubic');
                % Kriging model for prediction
                dmodel=dacefit(trainset_S, trainset_Y, 'regpoly1', 'corrgauss',...
                	ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                yPred_all(id_a:id_e,2)=predictor(validation_S,dmodel);
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end 
    end
%Mixture model: cubic RBF & MARS    
elseif strcmp(Surrogate,'MIX_RcM')  
    if ~forcekfold  % leave-1-out cross-validation
        for jj=1:m
            %cubic RBF model prediction
            [lambda, gamma]=RBF([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cubic'); 
            yPred_all(jj,1)=RBF_eval(Data.S(jj,:),[Data.S(1:jj-1,:);...
                Data.S(jj+1:end,:)],lambda,gamma,'cubic');
            % MARS model prediction
            mmodel = aresbuild([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)]); %MARS
            yPred_all(jj,2)=arespredict(mmodel,Data.S(jj,:)); 
        end
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);
        for jj=1:ceil(m/kfold)
            %determine validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %lef-tovers
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            %determine training set
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
                %RBF model prediction
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                yPred_all((jj-1)*kfold+1:jj*kfold,1)=RBF_eval(validation_S,...
                    trainset_S,lambda,gamma,'cubic');
                %MARS model prediction
                mmodel = aresbuild(trainset_S,trainset_Y); 
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=arespredict(mmodel,validation_S); 
            else %left-overs
                id_a=(jj-1)*kfold+1;
                if mod(m,kfold)>0
                    id_e=(jj-1)*kfold+mod(m,kfold);
                else
                    id_e=jj*kfold;
                end
                %RBF model prediction
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                yPred_all(id_a:id_e,1)=RBF_eval(validation_S,trainset_S,...
                    lambda,gamma,'cubic');
                %MARS model prediction
                mmodel = aresbuild(trainset_S,trainset_Y); 
                yPred_all(id_a:id_e,2)=arespredict(mmodel,validation_S); 
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end 
    end
%Mixture model: cubic RBF & full cuboc polynomial    
elseif strcmp(Surrogate,'MIX_RcPc')  
    if ~forcekfold  % leave-one-out cross-validation
        for jj=1:m
            %RBF model predictions
            [lambda, gamma]=RBF([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cubic'); 
            yPred_all(jj,1)=RBF_eval(Data.S(jj,:),[Data.S(1:jj-1,:);...
                Data.S(jj+1:end,:)],lambda,gamma,'cubic');
            % Polynomial model predictions
            beta=POLY([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cub'); 
            yPred_all(jj,2)=POLY_eval(Data.S(jj,:),beta,'cub');
        end
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);
        for jj=1:ceil(m/kfold)
            %validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            %training set
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
                % cubic polynomial model predictions
                beta=POLY(trainset_S,trainset_Y,'cub'); %Polynomial
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=POLY_eval(validation_S,beta,'cub');
            else %left-overs
                id_a=(jj-1)*kfold+1;
                if mod(m,kfold)>0
                    id_e=(jj-1)*kfold+mod(m,kfold);
                else
                    id_e=jj*kfold;
                end
                %RBF model prediction
                [lambda, gamma]=RBF(trainset_S,trainset_Y,'cubic');
                yPred_all(id_a:id_e,1)=RBF_eval(validation_S,trainset_S,...
                    lambda,gamma,'cubic');
                %cubic polynomial model prediction
                beta=POLY(trainset_S,trainset_Y,'cub'); %Polynomial
                yPred_all(id_a:id_e,2)=POLY_eval(validation_S,beta,'cub');
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end 
    end
%Mixture model: Kriging with Gaussian regression model& MARS    
elseif strcmp(Surrogate,'MIX_KgM')  
    if ~forcekfold  %leave-one-out cross-validation
        for jj=1:m
            %Kriging model predictions
            dmodel = dacefit([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)], 'regpoly1', ...
                'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
            if ~isstruct(dmodel)
                disp('Building the kriging model failed')
                return
            else
                yPred_all(jj,1)=predictor(Data.S(jj,:), dmodel);  
            end
            % MARS model predictions
            mmodel = aresbuild([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)]); 
            yPred_all(jj,2)=arespredict(mmodel,Data.S(jj,:)); 
        end    
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);
        for jj=1:ceil(m/kfold)
            %validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            %training set
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
                %Kriging model predictions
                dmodel = dacefit(trainset_S,trainset_Y, 'regpoly1', ...
                    'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                yPred_all((jj-1)*kfold+1:jj*kfold,1)=predictor(validation_S,dmodel);
                %MARS model predictions
                mmodel = aresbuild(trainset_S,trainset_Y); 
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=arespredict(mmodel,validation_S); 
            else %left-overs
                id_a=(jj-1)*kfold+1;
                if mod(m,kfold)>0
                    id_e=(jj-1)*kfold+mod(m,kfold);
                else
                    id_e=jj*kfold;
                end
                %Kriging model prediction
                dmodel = dacefit(trainset_S,trainset_Y, 'regpoly1', ...
                    'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                yPred_all(id_a:id_e,1)=predictor(validation_S,dmodel);
                %MARS model prediction
                mmodel = aresbuild(trainset_S,trainset_Y); 
                yPred_all(id_a:id_e,2)=arespredict(mmodel,validation_S); 
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end     
    end
%Mixture model: Kriging with Gaussian correlation model & full cubic
%polynomial
elseif strcmp(Surrogate,'MIX_KgPc')  
    if ~forcekfold  %leave-one-out cross-validation
        for jj=1:m
            %Kriging model prediction
            dmodel = dacefit([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)], 'regpoly1', ...
                'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
            if ~isstruct(dmodel)
                disp('Building the kriging model failed')
                return
            else
                yPred_all(jj,1)=predictor(Data.S(jj,:), dmodel);  
            end
            %full cubic polynomial model prediction
            beta=POLY([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cub'); 
            yPred_all(jj,2)=POLY_eval(Data.S(jj,:),beta,'cub');
        end    
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);   
        for jj=1:ceil(m/kfold)
            %validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            %validation set
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
                %Kriging model prediction
                dmodel = dacefit(trainset_S,trainset_Y, 'regpoly1', ...
                    'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                yPred_all((jj-1)*kfold+1:jj*kfold,1)=predictor(validation_S,dmodel);
                % Polynomial model prediction
                beta=POLY(trainset_S,trainset_Y,'cub');
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=POLY_eval(validation_S,beta,'cub');
            else %left-overs
                id_a=(jj-1)*kfold+1;
                if mod(m,kfold)>0
                    id_e=(jj-1)*kfold+mod(m,kfold);
                else
                    id_e=jj*kfold;
                end
                %Kriging model prediction
                dmodel = dacefit(trainset_S,trainset_Y, 'regpoly1', ...
                    'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                yPred_all(id_a:id_e,1)=predictor(validation_S,dmodel);
                % Polynomial model prediction
                beta=POLY(trainset_S,trainset_Y,'cub'); 
                yPred_all(id_a:id_e,2)=POLY_eval(validation_S,beta,'cub');
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end     
    end        
%Mixture model: Kriging with Gaussian correlation model & reduced cubic
%polynomial
elseif strcmp(Surrogate,'MIX_KgPcr')  
    if ~forcekfold  %leave-one-out cross-validation
        for jj=1:m
            %Kriging model prediction
            dmodel = dacefit([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)], 'regpoly1', ...
                'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
            if ~isstruct(dmodel)
                disp('Building the kriging model failed')
                return
            else
                yPred_all(jj,1)=predictor(Data.S(jj,:), dmodel);  
            end
            %full cubic polynomial model prediction
            beta=POLY([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'cubr'); 
            yPred_all(jj,2)=POLY_eval(Data.S(jj,:),beta,'cubr');
        end    
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);   
        for jj=1:ceil(m/kfold)
            %validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            %validation set
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
                %Kriging model prediction
                dmodel = dacefit(trainset_S,trainset_Y, 'regpoly1', ...
                    'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                yPred_all((jj-1)*kfold+1:jj*kfold,1)=predictor(validation_S,dmodel);
                % Polynomial model prediction
                beta=POLY(trainset_S,trainset_Y,'cubr'); 
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=POLY_eval(validation_S,beta,'cubr');
            else %left-overs
                id_a=(jj-1)*kfold+1;
                if mod(m,kfold)>0
                    id_e=(jj-1)*kfold+mod(m,kfold);
                else
                    id_e=jj*kfold;
                end
                %Kriging model prediction
                dmodel = dacefit(trainset_S,trainset_Y, 'regpoly1', ...
                    'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                yPred_all(id_a:id_e,1)=predictor(validation_S,dmodel);
                % Polynomial model prediction
                beta=POLY(trainset_S,trainset_Y,'cubr'); 
                yPred_all(id_a:id_e,2)=POLY_eval(validation_S,beta,'cubr');
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end     
    end        
%Mixture model: Kriging with Gaussian correlation model & reduced quadratic
%polynomial
elseif strcmp(Surrogate,'MIX_KgPqr')  
    if ~forcekfold  %leave-one-out cross-validation
        for jj=1:m
            %Kriging model prediction
            dmodel = dacefit([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)], 'regpoly1', ...
                'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));  
            if ~isstruct(dmodel)
                disp('Building the kriging model failed')
                return
            else
                yPred_all(jj,1)=predictor(Data.S(jj,:), dmodel);  
            end
            %full cubic polynomial model prediction
            beta=POLY([Data.S(1:jj-1,:);Data.S(jj+1:end,:)],...
                [Data.Ymed(1:jj-1,:);Data.Ymed(jj+1:end,:)],'quadr'); 
            yPred_all(jj,2)=POLY_eval(Data.S(jj,:),beta,'quadr');
        end    
    else %k-fold cross-validation
        Ss=Data.S(mix,:); %resorted sample sites
        Ys=Data.Ymed(mix);   
        for jj=1:ceil(m/kfold)
            %validation set
            if jj*kfold<=m
                validation_S=Ss((jj-1)*kfold+1:jj*kfold,:);
            else %left-overs
                validation_S=Ss((jj-1)*kfold+1:end,:);
            end
            %validation set
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
                %Kriging model prediction
                dmodel = dacefit(trainset_S,trainset_Y, 'regpoly1', ...
                    'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                yPred_all((jj-1)*kfold+1:jj*kfold,1)=predictor(validation_S,dmodel);
                % Polynomial model prediction
                beta=POLY(trainset_S,trainset_Y,'quadr'); %Polynomial
                yPred_all((jj-1)*kfold+1:jj*kfold,2)=POLY_eval(validation_S,beta,'quadr');
            else %left-overs
                id_a=(jj-1)*kfold+1;
                if mod(m,kfold)>0
                    id_e=(jj-1)*kfold+mod(m,kfold);
                else
                    id_e=jj*kfold;
                end
                %Kriging model prediction
                dmodel = dacefit(trainset_S,trainset_Y, 'regpoly1', ...
                    'corrgauss', ones(1,Data.dim), 1e-1*ones(1,Data.dim), 20*ones(1,Data.dim));
                yPred_all(id_a:id_e,1)=predictor(validation_S,dmodel);
                % Polynomial model prediction
                beta=POLY(trainset_S,trainset_Y,'quadr'); 
                yPred_all(id_a:id_e,2)=POLY_eval(validation_S,beta,'quadr');
                yPred_all=sortrows([yPred_all,mix(:)],Data.numberOfModels+1);
                yPred_all=yPred_all(:,1:Data.numberOfModels);
            end
        end     
    end        
        
%add other 2-model mixtures here          
end
        
%% characteristics of single models 
Corcoef=cc_calc(Data.Ymed,yPred_all); %correlation coefficients
RMSE=RMSE_calc(Data.Ymed,yPred_all); %root mean squared errors
MAE=MAE_cal(Data.Ymed,yPred_all); %maximum absolute error
MAD=MAD_cal(Data.Ymed,yPred_all); %median absolute deviation
   
%Scale correlation coefficients
Corcoef_t=Corcoef;
if any(Corcoef_t<0) %scale correlation coefficients to interval [0,1]
    Corcoef_t=(Corcoef_t-min(Corcoef_t))./max(Corcoef_t-min(Corcoef_t));
end
if all(isnan(Corcoef)) %all models get same CC
    Corcoef_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    Corcoef_t(Corcoef_t==0)=0.000001; %every CC most be positive
    Corcoef_t=Corcoef_t(:)./(sum(Corcoef_t(:)));
end
%scale root mean squared errors
if all(RMSE==0) %Same RMSE for each model
    RMSE_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    RMSE_t=(1./RMSE(:))./sum(1./RMSE(:));
end
%scale MAximum absolute errors
if all(MAE==0) %Same MAE for each model
    MAE_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    MAE_t=(1./MAE(:))./sum(1./MAE(:));
end
%scale median absolute deviation
if all(MAD==0) %Same MAD for each model
    MAD_t=1/Data.numberOfModels*ones(Data.numberOfModels,1);
else
    MAD_t=(1./MAD(:))./(sum(1./MAD));
end

%compute basic probability assignments (BPA)    
comb=[1 NaN; NaN 2]; %individual models
OVpcr=[comb Corcoef_t(:), RMSE_t(:), MAE_t(:), MAD_t(:)]; %evidences
    
intround=true; 
%DST
for yy =2:Data.numberOfModels
    if intround %first model
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
%compute prediction of mixture model
Ymixtest=zeros(size(yPred_all,1),nchoosek(Data.numberOfModels,2));
for ii = 1:nchoosek(Data.numberOfModels,2)
    Ymixtest(:,ii)=Intervals(1,Data.numberOfModels+1)*yPred_all(:,1)+...
        Intervals(1,Data.numberOfModels+1)*yPred_all(:,2); %mixture model predictions at sample points 
end
%compute model characteristics for mixture model
%compute correlation coefficient of prediction by mixture    
Corcoef=[Corcoef;cc_calc(Data.Ymed,Ymixtest)];
if any(Corcoef<0)
    Corcoef=(Corcoef-min(Corcoef))./max(Corcoef-min(Corcoef));
end
Corcoef(Corcoef==0)=0.000001;
RMSE=[RMSE;RMSE_calc(Data.Ymed,Ymixtest)]; %compute RMSE of mixture
MAE=[MAE;MAE_cal(Data.Ymed,Ymixtest)]; %compute MAE of mixture
MAD=[MAD;MAD_cal(Data.Ymed,Ymixtest)]; % compute MAD of mixture     


%% Combine pure models to mixed model     
%all available model combinations (pure models and 2-model mixture)
comb=[1 NaN; NaN 2; 1 2]; 

%matrix containing model combination and corresponding evidences (high numbers mean high evidence)
OVpcr=[comb Corcoef(:)./(sum(Corcoef(:))) (1./RMSE(:))./sum(1./RMSE(:))...
    (1./MAE(:))./sum(1./MAE(:)) (1./MAD(:))./(sum(1./MAD))];
clear Corcoef; clear RMSE; clear MAE; clear MAD
clear Corcoef_t; clear RMSE_t; clear MAE_t; clear MAD_t;

intround=true;
%DST to determine weights in models
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
w=Intervals([2:3],end);%weights_combi(Intervals,ModelParam.newmodel);
modelweights=w/sum(w); %scale model weights so that w1+w2=1

end % function