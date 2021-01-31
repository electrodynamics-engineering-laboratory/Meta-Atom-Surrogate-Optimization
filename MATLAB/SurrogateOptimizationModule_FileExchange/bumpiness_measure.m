function hn=bumpiness_measure(x,Data,flag, target, tolerance, lambda, gamma)

%function that computes the "bumpiness" of the response surface given the
%currnt objective function value target.
%--------------------------------------------------------------------------
%Copyright (c) 2012 by Juliane Mueller
%
% This file is part of the surrogate model module toolbox.
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%Tampere University of Technology, Finland
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%input:
%x - optimization variable, the point that minimizes the bumpiness of the
%    surface
%Data - structure array containing all problem information
%flag - string determining which radial basis function to use
%target - target for objective function value 
%tolerance - scalar determining when two sample sites are too
%            close/coincide
%lambda, gamma - parameter values for the response surface given all sample
%                data obtained so far
%
%Output: 
%hn - bumpiness measure
%--------------------------------------------------------------------------

m=size(Data.S,1); %dimension of sample site matrix
P=[Data.S, ones(m,1)];%establish matrix of sample sites
R_y=sqrt(sum((Data.S-repmat(x,m,1)).^2,2)); %compute distance between x and all already sampled points

if any(R_y < tolerance) %point x is too close to already sampled points
    hn=0; %give the bumpiness a bad bumpiness function value 
else
    %pairwise distances between sample points
    R=zeros(m,m); %initialize distance matrix
    for ii = 1:m 
        for jj = ii:m
            R(ii,jj)=sum((Data.S(ii,:)-Data.S(jj,:)).^2,2);
            R(jj,ii)=R(ii,jj);
        end
    end
    R=sqrt(R); %R contains the pairwise euclidean distances between all points
    %compute radial basis function values
    if strcmp(flag,'cubic') %cubic RBF
        Phi= R.^3;
        Phi_y=R_y.^3; %value for the point x
        m0=1;
    elseif strcmp(flag,'TPS') %thin plate spline RBF
        R(R==0)=1;
        Phi=R.^2.*log(R);
        Phi_y=R_y.^2.*log(R_y); %value for the point x
        m0=1;
    elseif strcmp(flag,'linear') %linear RBF
        Phi=R;
        Phi_y=R_y; %value for the point x
        m0=0;
    end
    %build augmented matrices by adding the point x to the current matrices
    Phi_aug=[Phi,Phi_y;Phi_y', 0]; %augmented radial basis function value matrix
    P_aug=[P;[x, 1]]; %augmented sample site matrix
    A=[Phi_aug, P_aug; P_aug', zeros(Data.dim+1,Data.dim+1)]; %matrix A*v = rhs (see Gutmann 2001)
    rhs=[zeros(m,1); 1; zeros(Data.dim+1,1)]; %(m+1)th variable corresponds to mu
    v=A\rhs; %this one requires O(n^3) operations, can be modified to require only O(n^2)
    mu = v(m+1); %pick the value for mu from the solution
    %mu must be positive for TPS and cubic RBF, and negative for linear RBF
    if abs(mu) <1e-6
        mu=0;
    elseif (mu<0 && (strcmp(flag,'TPS') || strcmp(flag,'cubic'))) || ...
            (mu>0 && (strcmp(flag,'linear')))
        error('mu in bumpiness_measure.m')
    end
    yhat=RBF_eval(x,Data.S,lambda,gamma,flag); %predicted function value at site x
    gn=(-1)^(m0+1)*mu*(yhat-target)^2; %bumpiness measure
    hn=-1/gn; %minimize -1/gn to avoid numerical difficulties at already sampled points
end

end%function