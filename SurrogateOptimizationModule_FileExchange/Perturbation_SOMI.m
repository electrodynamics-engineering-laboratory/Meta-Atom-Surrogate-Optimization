function CandPoint = Perturbation_SOMI(Data,NCandidates,stdev_int,stdev_cont,P)

%              
%function that creates candidate points for the next sample site. 
%four groups are generated: 
%a) only continuous perturbation (small, medium, large)
%b) only integer perturbation (small, medium, large)
%c) continuous and integer perturbation (small, medium, large)
%d) randomly generated points over whole domain
     
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
%Data - structure, containing all problem information
%NCandidates - integer, number of candidate points in each group
%stdev_int - vector, contains perturbation ranges for integers
%stdev_cont - vector, contains perturbation ranges for continuous variables
%P - scalar, perturbation probability
%--------------------------------------------------------------------------

% sequential execution:
%Group1: perturb only continuous variables  
d_cont=length(Data.continuous);   %dimension of continuous variables

%set candidate points same as xbest, and then perturb only continuous
%variables
CandPoint_G1=repmat(Data.xbest,NCandidates,1);
for ii = 1:NCandidates %for each canidate point
    for jj=1:d_cont %for each continuous variables
        if rand(1)<= P %check if perturbation wanted
            r=randperm(length(stdev_cont));
            CandPoint_G1(ii,Data.continuous(jj)) = max(Data.xlow(Data.continuous(jj)),...
                min(CandPoint_G1(ii,Data.continuous(jj))+stdev_cont(r(1))*randn(1),...
                Data.xup(Data.continuous(jj))));
        end
    end
end


%Group2: perturb only integer variables     
%set candidate points same as xbest, and then perturb only integer
%variables
CandPoint_G2=repmat(Data.xbest,NCandidates,1);
for ii = 1:NCandidates %for each canidate point
    for jj=1:Data.dim %for each integer variables
        if ismember(jj,Data.integer) && rand(1)<= P %check if perturbation wanted
            r=randperm(length(stdev_int));  %select random perturbation range
            p=randn(1);
            sign_p=sign(p);
            if sign_p <0
                perturbation= min(-1,round(stdev_int(r(1))*randn(1)));
            else
                perturbation= max(1,round(stdev_int(r(1))*randn(1)));
            end
            CandPoint_G2(ii,jj) = max(Data.xlow(jj),...
                min(CandPoint_G2(ii,jj)+perturbation,...
                Data.xup(jj)));
        end
    end
end

%Group 3
%set candidate points same as xbest, and then perturb integer & continuous
%variables
CandPoint_G3=repmat(Data.xbest,NCandidates,1);
for ii = 1:NCandidates %for each canidate point
    for jj=1:Data.dim %for each variables
        if ismember(jj,Data.integer) 
            if rand(1)<= P %check if perturbation wanted
                r=randperm(length(stdev_int));
                p=randn(1);
                sign_p=sign(p);
                if sign_p <0
                    perturbation= min(-1,round(stdev_int(r(1))*randn(1)));
                else
                    perturbation= max(1,round(stdev_int(r(1))*randn(1)));
                end
                CandPoint_G3(ii,jj) = max(Data.xlow(jj),...
                    min(CandPoint_G3(ii,jj)+perturbation,Data.xup(jj)));
            end
        else
            if rand(1)<= P %check if perturbation wanted
                r=randperm(length(stdev_cont));
                CandPoint_G3(ii,jj) = max(Data.xlow(jj),...
                    min(CandPoint_G3(ii,jj)+stdev_cont(r(1))*randn(1),...
                    Data.xup(jj)));
            end
        end
    end
end


%Group4: uniformly sampled points over variable domain
CandPoint_G4=repmat(Data.xlow,NCandidates,1) + rand(NCandidates,length(Data.xlow)).*...
    repmat(Data.xup-Data.xlow,NCandidates,1);
%round integer variables
CandPoint_G4(:,Data.integer)=round(CandPoint_G4(:,Data.integer));

%put all candidate points and predicted objective function values together
CandPoint=[CandPoint_G1;CandPoint_G2;CandPoint_G3; CandPoint_G4];

%delete unnecessary data
clear CandPoint_G1 CandPoint_G2 CandPoint_G3 CandPoint_G4 ...

end %function