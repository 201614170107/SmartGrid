function [Popt,Qopt,rhopt]=ComputePdist(mpc)
VS              =   ComputeLoadLess(mpc);
[Matsx,Matsy]   =   MakeMats(mpc);
nb              =   size(mpc.bus,1);
nsb             =   find(~(mpc.bus(:,2)==3));
sb              =   find((mpc.bus(:,2)==3));
n               =   size(Matsx,1);
k               =   size(VS,2);
rhos            =   [0,1];
Mo              =   sum(Matsx,3)*real(VS(sb,1))+sum(Matsy,3)*imag(VS(sb,1));
Popt            =   [];
Qopt            =   [];
rhopt           =   0;

for kt=1:6
    rho             =   mean(rhos);
    [c,P,Q]         =   OptPdist(ComputeMs(rho),Mo,0);
    if(c>0.1)
        rhos    =   [rho,rhos(2)];
        rhopt   =   rho;
        Popt    =   P;
        Qopt    =   Q;
    else
        rhos    =   [rhos(1),rho];
    end
end

    function [Ms]=ComputeMs(rho)
        Ms  =   zeros(n,n,k);
        for it=1:k
            x=real(rho*VS(:,it)-rho*VS(sb,1));
            y=imag(rho*VS(:,it)-rho*VS(sb,1));
            Ms(:,:,it)=sum(bsxfun(@times,Matsx,permute(x,[3,2,1])),3)+...
                sum(bsxfun(@times,Matsy,permute(y,[3,2,1])),3);
        end                
    end

Sols    =   [rhopt*VS+(1-rhopt)*VS(sb,1),VS(sb,1)*ones(nb,1)];

Sols    =   [real(Sols(nsb,:));imag(Sols(nsb,:))];

end