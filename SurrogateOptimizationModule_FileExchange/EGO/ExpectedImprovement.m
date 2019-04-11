function ExpImp = ExpectedImprovement(Data,x)

%EI calculates expected improvement for EGO algorithm
%implementation based on code by Forrester
%
%Input:
%Data - structure,  contains information about the problem
%x - point at which expected improvement is computed
%
%Output:
%ExpImp - expected improvement at the point x


yBest=min(Data.Ymed); %use median-corrected data
theta=10.^Data.Theta;

m=size(Data.S,1); % number of sample points
e=ones(m,1);      % vector of ones

try
    mu=(e'*(Data.U\(Data.U'\Data.Ymed)))/(e'*(Data.U\(Data.U'\e))); %calculate mu

    ssq=((Data.Ymed-e*mu)'*(Data.U\(Data.U'\(Data.Ymed-e*mu))))/m; % calculate sigma^2
    f=zeros(size(x,1),1);
    Ssqr=zeros(size(f));
    ExpImp=zeros(size(f));

    for kk=1:size(x,1)
        psi=exp(-sum(repmat(theta,size(Data.S,1),1).*abs(Data.S-repmat(x(kk,:),size(Data.S,1),1)).^2,2)); %correlation model
        f(kk,1)=mu+psi'*(Data.U\(Data.U'\(Data.Ymed-e*mu)));   %calculate predictions
        Ssqr(kk,1)=ssq*(1-psi'*(Data.U\(Data.U'\psi)));     %get error


        if Ssqr(kk)==0
            ExpImp(kk,1)=0;
        else
            ei_termone=(yBest-f(kk))*(0.5+0.5*erf((1/sqrt(2))*((yBest-f(kk))/sqrt(abs(Ssqr(kk))))));
            ei_termtwo=sqrt(abs(Ssqr(kk)))*(1/sqrt(2*pi))*exp(-(1/2)*((yBest-f(kk))^2/Ssqr(kk)));
            ExpImp(kk,1)=ei_termone+ei_termtwo;
        end

    end
catch %in case matrix is ill-conditioned, Data.U will not be of correct size
    disp('correlation matrix not symm pos def')
       mu=e'*(Data.Psi\Data.Ymed)/(e'*(Data.Psi\e));
       ssq=(Data.Ymed-e*mu)'*(Data.Psi\(Data.Ymed-e*mu))/m;
       f=zeros(size(x,1),1);
       Ssqr=zeros(size(f));
       ExpImp=zeros(size(f));
       for kk=1:size(x,1)
            psi=exp(-sum(repmat(theta,size(Data.S,1),1).*abs(Data.S-repmat(x(kk,:),size(Data.S,1),1)).^2,2));
            f(kk,1)=mu+psi'*(Data.Psi\(Data.Ymed-e*mu));
            Ssqr(kk,1)=ssq*(1-psi'*(Data.Psi\psi)); 
            %Ssqr(kk,1)=ssq*(1-psi'*(Data.Psi\psi)+(1-e'*(Data.Psi\psi))/(e'*(Data.Psi\e)));
            if Ssqr(kk)==0
                ExpImp(kk,1)=0;
            else
                ei_termone=(yBest-f(kk))*(0.5+0.5*erf((1/sqrt(2))*((yBest-f(kk))/sqrt(abs(Ssqr(kk))))));
                ei_termtwo=sqrt(abs(Ssqr(kk)))*(1/sqrt(2*pi))*exp(-(1/2)*((yBest-f(kk))^2/Ssqr(kk)));
                ExpImp(kk,1)=ei_termone+ei_termtwo;
            end
       end
end    
end%function