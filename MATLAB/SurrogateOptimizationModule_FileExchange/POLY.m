function b=POLY(S,Y,flag)
%computes the coefficients of a regression polynomial of the specified type
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
%S - sample site matrix
%Y - objective function values corresponding to the points in S
%flag - determines the order of polynomial model
%
%Output:
%b - vector with parameters
%--------------------------------------------------------------------------

[m,n]=size(S);
%linear
if strcmp(flag,'lin') %first order polynomial
    X=[ones(m,1) S];
    b=((X'*X)\X')*Y;  
%full quadratic    
elseif strcmp(flag,'quad') % f= b0 + b1 x1 + b2 x2+...+bn  xn + b12 x1x2 + b13 x1x3+...+b11 x1^2 + ...+bnn xn^2 
    X=[ones(m,1) S S.^2];
    ii=1;
    while ii < n
        jj=ii+1;
        while jj <= n
            x=S(:,ii).*S(:,jj);
            jj=jj+1;
            X=[X,x];
        end
        ii=ii+1;
    end
    b=((X'*X)\X')*Y; 
%full cubic    
elseif strcmp(flag,'cub') % f= b0 + b1x1 + b1b2 x1x2 + b1b2b3 x1x2x3+...+b11 x1^2 + ...+bnnn xn^3 
    X=[ones(m,1) S S.^2 S.^3];
    ii=1;
    while ii < n %x_i x_j type terms
        jj=ii+1;
        while jj <= n            
            x=S(:,ii).*S(:,jj);
            jj=jj+1;
            X=[X,x];
        end
        ii=ii+1;
    end
    
    ii=1;
    while ii < n-1 %x_i x_j x_k type terms
        jj=ii+1;
        while jj < n
            kk= jj +1; 
            while kk <= n
                x=S(:,ii).*S(:,jj).*S(:,kk);
                kk=kk+1;
                X=[X,x];
            end
            jj=jj+1;
        end
        ii=ii+1;
    end   
    b=((X'*X)\X')*Y;     
%reduced quadratic    
elseif strcmp(flag,'quadr') %f = b0+ b1 x1 + b2 x2 + ...  bn xn + b11 x1^2 + b22 x2^2...
    X=[ones(m,1) S S.^2];
    b=((X'*X)\X')*Y;     
%reduced cubic    
elseif strcmp(flag,'cubr') %f = b0+ b1 x1 + b2 x2 + ...  bn xn + b11 x1^2 + b22 x2^2...+b111 x1^3 + b222 x2^3+...
    X=[ones(m,1) S S.^2 S.^3];
    b=((X'*X)\X')*Y;     
end

end %function