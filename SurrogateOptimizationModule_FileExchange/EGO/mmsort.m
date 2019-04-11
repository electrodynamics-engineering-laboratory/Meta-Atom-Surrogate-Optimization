function Index = mmsort(X3D,p)
% Ranks sampling plans according to the Morris-Mitchell criterion
% definition. Note: similar to phisort, which uses the numerical quality
% criterion Phiq as a basis for the ranking.
%
% Inputs:
%       X3D - three-dimensional array containing the sampling plans to be
%       ranked.
%       p - the distance metric to be used (p=1 rectangular - default, p=2
%       Euclidean)
%
% Output:
%       Index - index array containing the ranking
%
% Copyright 2007 A Sobester
%

if ~exist('p','var')
    p = 1;
end

% Pre-allocate memory
Index = (1:size(X3D,3));

% Bubble-sort
swap_flag = 1;

while swap_flag==1
    swap_flag = 0;
    i=1;
    while i<=length(Index)-1
        if mm(X3D(:,:,Index(i)),X3D(:,:,Index(i+1)),p)==2
            buffer = Index(i);
            Index(i) = Index(i+1);
            Index(i+1) = buffer;
            swap_flag = 1;
        end
        i = i + 1;
    end
end
