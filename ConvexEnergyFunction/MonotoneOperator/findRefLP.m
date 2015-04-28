function [vref,Popt,Pm,rhopt]=findRefLP(mpc,varargin)
[Matsx,Matsy]=MakeMats(mpc);
n=size(Matsx,1);
m=size(Matsx,3);
nsb=~(mpc.bus(:,2)==3);
sb=find(mpc.bus(:,2)==3);

Vin=ones(m,1);
if(nargin>1)
    Vin=varargin{1};
end


th=0;%linspace(0,2*pi,1);
thopt=0;
Popt=[];

Ms=cat(3,Matsx,Matsy);
Pm=amd(sum(Ms,3));
Ms=Ms(Pm,Pm,:);

rhopt=0;

nb=size(mpc.branch,1);
Amat=sparse(mpc.branch(:,1),1:nb,1,n,nb)-sparse(mpc.branch(:,2),1:nb,1,n,nb);


for it=1:length(th)
    V   =   Vin*exp(1i*th(it));
    x   =   real(V);
    y   =   imag(V);
    M   =   sum(Ms,3)*x(1)+sum(Ms,3)*y(1);
    Mn  =   Ms(:,:,[nsb;nsb]);
    A   =   Amat(nsb,:);
    if(rhopt>0)
    chk =   CheckRhoLP(rhopt,M,Mn,A);
    if(~chk)
        continue;
    end    
    end
    
    rhos=   [rhopt,3];
    for kt=1:5
        chk =   CheckRhoLP(mean(rhos),M,Mn,A);
        if(chk)
            rhos=[mean(rhos),max(rhos)];
        else
            rhos=[min(rhos),mean(rhos)];
        end
    end
    
    rho=min(rhos);
    if(rho>rhopt)
        rhopt=rho;
         [~,P]=CheckRhoLP(rhopt,M,Mn,A); 
        Popt=P*sparse(1:n,Pm,1,n,n);
        thopt=th(it);
    end
    disp(rhopt);
end

  

vref=[cos(thopt);sin(thopt)];

end