function samplesites=cornerpoints(Data,number_startpoints)

%generates initial experimental design with points in corners and center
%point
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
%Data - structure with problem information
%number_startpoints - number of desired points in starting design
%
%output:
%samplesites - matrix with points in initial experimental design
%--------------------------------------------------------------------------

if number_startpoints > 2^Data.dim+1
    error('There are not enough corners to create the requested number of starting points')
end

S=zeros(2^Data.dim,Data.dim); %initialize sample site matrix 
for ii = 1:Data.dim
    S(:,ii)=repmat([repmat(Data.xlow(ii),2^(Data.dim-ii),1);...
        repmat(Data.xup(ii),2^(Data.dim-ii),1)],2^(ii-1),1);
end


%if there are more corner points than desired number of starting  points,
%randomly select corner points
if size(S,1) >number_startpoints-1 
    r=randperm(size(S,1));
    use=S(r(1:number_startpoints-1),:); %select randomly corners
    samplesites=[use;Data.xlow+0.5*(Data.xup-Data.xlow)]; %add center point
else
    samplesites=[S;Data.xlow+0.5*(Data.xup-Data.xlow)];%add center point
end

end %function