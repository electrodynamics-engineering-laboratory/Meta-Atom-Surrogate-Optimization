function [NegLnLike,Psi,U]=likelihood_new(Data,x)
%
% Calculates the negative of the concentrated ln-likelihood
%
% Inputs:
%x - vector of log(theta) parameters
%Data - structure containing information about the optimization probelms
%
%Outputs:
%NegLnLike - concentrated log-likelihood *-1 for minimising
%Psi - correlation matrix
%U - Choleski factorisation of correlation matrix
%
% Copyright 2007 A I J Forrester
%

theta=10.^x;
n=size(Data.S,1);
one=ones(n,1);

D=repmat(theta,size(Data.dist,1),1).*Data.dist;
dist=exp(-sum(D,2));

dd=sortrows([Data.ij,dist],2);

Psi=triu(ones(n),1);
Psi(Psi==1)=dd(:,end);
Psi=Psi+Psi'+eye(n)+eye(n).*eps;


% Cholesky factorisation
[U,p]=chol(Psi);

% Use penalty if ill-conditioned
if p>0
    NegLnLike=1e4;
else
    
    % Sum lns of diagonal to find ln(abs(det(Psi)))
    LnDetPsi=2*sum(log(abs(diag(U))));

    % Use back-substitution of Cholesky instead of inverse
    mu=(one'*(U\(U'\Data.Ymed)))/(one'*(U\(U'\one)));
    SigmaSqr=((Data.Ymed-one*mu)'*(U\(U'\(Data.Ymed-one*mu))))/n;
    NegLnLike=-1*(-(n/2)*log(SigmaSqr) - 0.5*LnDetPsi);
end
end %function