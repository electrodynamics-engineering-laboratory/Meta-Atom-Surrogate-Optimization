function Yest=POLY_eval(X,b,flag)

%uses polynomial regression model to predict objective function values
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
%Input 
%S - matrix containing points where prediction wanted
%b - parameter vector computed with POLY.m
%flag - order of polynomial. must be same as flag  used for computing vector b
%
%Output
%Yest - predicted objective function values corresponding to the points in
%S
%--------------------------------------------------------------------------

[m,n]=size(X);

if strcmp(flag,'lin') %first order polynomial
    Xs=[ones(m,1) X];
    Yest=Xs*b;    
elseif strcmp(flag,'quad') % f= b0 + b1 x1 + b2 x2+...+bn  xn + b12 x1x2 + b13 x1x3+...+b11 x1^2 + ...+bnn xn^2 
    Xs=[ones(m,1) X X.^2];
    ii=1;
    while ii < n
        jj=ii+1;
        while jj <= n
            x=X(:,ii).*X(:,jj);
            jj=jj+1;
            Xs=[Xs,x];
        end
        ii=ii+1;
    end
    Yest=Xs*b;
elseif strcmp(flag,'cub') % f= b0 + b1x1 + b1b2 x1x2 + b1b2b3 x1x2x3+...+b11 x1^2 + ...+bnnn xn^3 
    Xs=[ones(m,1) X X.^2 X.^3];
    ii=1;
    while ii < n %x_i x_j type terms
        jj=ii+1;
        while jj <= n            
            x=X(:,ii).*X(:,jj);
            jj=jj+1;
            Xs=[Xs,x];
        end
        ii=ii+1;
    end
    
    ii=1;
    while ii < n-1 %x_i x_j x_k type terms
        jj=ii+1;
        while jj < n
            kk= jj +1; 
            while kk <= n
                x=X(:,ii).*X(:,jj).*X(:,kk);
                kk=kk+1;
                Xs=[Xs,x];
            end
            jj=jj+1;
        end
        ii=ii+1;
    end   
    Yest=Xs*b;
elseif strcmp(flag,'quadr') %f = b0+ b1 x1 + b2 x2 + ...  bn xn + b11 x1^2 + b22 x2^2...
    Xs=[ones(m,1) X X.^2];
    Yest=Xs*b;
elseif strcmp(flag,'cubr') %f = b0+ b1 x1 + b2 x2 + ...  bn xn + b11 x1^2 + b22 x2^2...+b111 x1^3 + b222 x2^3+...
    Xs=[ones(m,1) X X.^2 X.^3];
    Yest=Xs*b;
end
end %function