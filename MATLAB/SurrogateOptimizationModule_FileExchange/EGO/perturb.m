function X = perturb(X, PertNum)
% Interchanges pairs of randomly chosen elements within randomly chosen
% columns of a sampling plan a number of times. If the plan is a Latin
% hypercube, the result of this operation will also be a Latin hypercube.
%
% Inputs:
%       X - sampling plan
%       PertNum - the number of changes (perturbations) to be made to X.
%
% Output:
%       X - perturbed sampling plan
%
% Copyright 2007 A Sobester
%

if ~exist('PertNum','var')
	PertNum = 1;
end


[n,k] = size(X);

for  pert_count=1:PertNum
	col = floor(rand*k)+1;
	
	% Choosing two distinct random points
	el1 = 1; el2 = 1;	
	while el1==el2
		el1 = floor(rand*n)+1;
		el2 = floor(rand*n)+1;
	end

	%Swap the two chosen elements
	buffer = X(el1,col);
	X(el1,col) = X(el2,col);
	X(el2,col) = buffer;
end