function Phiq = mmphi(X, q, p)
% Calculates the sampling plan quality criterion of Morris and Mitchell.
%
% Inputs:
%       X - sampling plan
%       q - exponent used in the calculation of the metric
%       p - the distance metric to be used (p=1 rectangular - default, p=2
%           Euclidean)
%
% Output:
%       Phiq - sampling plan `space-fillingness' metric
%
% Copyright 2007 A Sobester
%

% Assume defaults if arguments list incomplete
if ~exist('p','var')
    p = 1;
end

if ~exist('q','var')
    q = 2;
end

% Calculate the distances between all pairs of
% points (using the p-norm) and build multipli-
% city array J
[J,d] = jd(X,p);

% The sampling plan quality criterion
Phiq = sum(J.*(d.^(-q)))^(1/q);