function intersections=dempster_rule(comp,OVnew,Data)
%dempster's rule of combination
%iteratively including next evidence set
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
%comp - matrix with models to combine and their evidences
%OVnew - matrix with model and their evidences from other criterion 
%Data - structure with all model information
%Output:
%intersections - "accumulated" evidences for each model in combination
%--------------------------------------------------------------------------

found=false;
intersections=[];
DR=1;

%compute conflict mass DR
for kk = 1:size(comp,1)%line 1
    for ll = 1:size(OVnew,1)%line 2
        if isempty(intersect(comp(kk,1:Data.numberOfModels),OVnew(ll,1:Data.numberOfModels))) &&...
                (~isnan(comp(kk,Data.numberOfModels+1))&&~isnan(OVnew(ll,Data.numberOfModels+1)))
            DR=DR-(comp(kk,Data.numberOfModels+1)*OVnew(ll,Data.numberOfModels+1));
        end
    end
end

b=(1 : Data.numberOfModels);

%compare comp to all other focal elements in OVnew and assign basic
%probabilities to corresponding sets
for pp=1:size(comp,1)
    for kk=1:size(OVnew,1)
        check=intersect(comp(pp,1:Data.numberOfModels),OVnew(kk,1:Data.numberOfModels)); %check if intersection is empty
        if ~isempty(check)
            addval=1/DR*(comp(pp,Data.numberOfModels+1)*OVnew(kk,Data.numberOfModels+1));
        else 
            addval=0;
        end
        I=NaN*ones(1,Data.numberOfModels);
        I(check)=b(check);
        if isempty(intersections)
            intersections=[I addval];
        elseif isempty(check)
            continue;
        else
            for ll =1:size(intersections,1)
                if ((sum(isnan(intersections(ll,1:Data.numberOfModels)))==sum(isnan(I))) &&...
                        (all(~isnan(intersections(ll,1:Data.numberOfModels))-~isnan(I)==0)) )|| all(isnan(I))
                    found=true;
                    if ~isnan(addval) 
                        intersections(ll,end)=intersections(ll,end)+addval;
                    end
                    break;
                end
            end
            if ~found
                intersections=[intersections; I addval];
            else
                found=false;
            end
        end
    end
end

%final intersections are basic probability assignments
%BPA=intersections;

end %function