function X = bestlh(n,k, Population, Iterations)
% Generates an optimized Latin hypercube by optimizing the Morris-Mitchell
% criterion for a range of exponents and plots the first two dimensions of
% the current hypercube throughout the optimization process.
%
% Inputs:
%       n - number of points required
%       k - number of design variables
%       Population - number of individuals in the evolutionary operation
%                    optimizer
%       Iterations - number of generations the evolutionary operation
%                    optimizer is run for
%       Note: high values for the two inputs above will ensure high quality
%       hypercubes, but the search will take longer.
%
% Output:
%       X - optimized Latin hypercube
%
% Copyright 2007 A Sobester
%

if k<2
	error('Latin hypercubes are not defined for k<2');
end

% List of qs to optimize Phi_q for
q = [1 2 5 10 20 50 100];

% Set the distance norm to rectangular for a faster search. This can be 
% changed to p=2 if the Euclidean norm is required.
p = 1;

% We start with a random Latin hypercube
XStart = rlh(n,k,1);

% For each q optimize Phi_q
for i=1:length(q)
    X3D(1:n,1:k,i) = mmlhs(XStart, Population, Iterations, q(i));
end

% Sort according to the Morris-Mitchell criterion
Index = mmsort(X3D,p);

% And the Latin hypercube with the best space-filling properties is...
X = X3D(:,:,Index(1));

end%function