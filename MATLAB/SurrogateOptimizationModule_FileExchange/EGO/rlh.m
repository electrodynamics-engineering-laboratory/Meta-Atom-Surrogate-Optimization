function X = rlh(n, k, Edges)
% Generates a random Latin hypercube within the [0,1]^k hypercube.
%
% Inputs:
%       n - desired number of points
%       k - number of design variables (dimensions)
%       Edges - if Edges=1 the extreme bins will have their centres on the
%               edges of the domain, otherwise the bins will be entirely 
%               contained within the domain (default setting). 
%
% Output:
%       X - Latin hypercube sampling plan of n points in k dimensions.
%
% Copyright 2007 A Sobester
%

if nargin < 3
    Edges = 0;
end

% Pre-allocate memory
X = zeros(n,k);

for i=1:k
   X(:,i) = randperm(n)';
end

if Edges == 1
    X = (X-1)/(n-1);
else
    X = (X-0.5)/n;
end