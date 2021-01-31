function [D,ij]=distanceupdate(Data)
%computes pairwise distances between sample sites
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
% Data - structure containing information about the problem
%
%Output:
%D - matrix with pairwise distances between variable values of points in S
%ij - matrix with indices of points for which pairwise distances computed
%--------------------------------------------------------------------------

p=2;  % added p definition (February 10)
n=size(Data.S,1);

% Pre-allocate memory

mzmax = n*(n-1) / 2;
ij = zeros(mzmax, 2);       % initialize matrix with indices
D = zeros(mzmax, Data.dim);        % initialize matrix with distances
ll = 0;
for k = 1 : n-1
  ll = ll(end) + (1 : n-k);
  ij(ll,:) = [repmat(k, n-k, 1) (k+1 : n)']; % indices for sparse matrix
  D(ll,:) = abs(repmat(Data.S(k,:), n-k, 1) - Data.S(k+1:n,:)); % differences between points
end

D=D.^p;

end %function