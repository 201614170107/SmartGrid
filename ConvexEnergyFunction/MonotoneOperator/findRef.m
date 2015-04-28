function [vref,Popt,Qopt,rhopt]=findRef(fname,opt,vf)
mpc=loadcase(fname);
fname=strcat(fname,'Popt.mat');
[Matsx,Matsy]=MakeMats(mpc);
n=size(Matsx,1);
m=size(Matsx,3);


nsb=~(mpc.bus(:,2)==3);

nb=size(mpc.branch,1);
Amat=sparse(mpc.branch(:,1),1:nb,1,n,nb)-sparse(mpc.branch(:,2),1:nb,1,n,nb);
A   =   Amat(nsb,:);
%Achk=	sum(bsxfun(@times,Matsx,permute(real(Vref),[3,2,1])),3)+...
 %           sum(bsxfun(@times,Matsy,permute(imag(Vref),[3,2,1])),3);
if(opt>0)
    th=0;
else
    th=linspace(0,2*pi,30);
end


    function [ch,Pop]=ffun(anyx,anyM,anyMs,anyNs)
        if(opt>0)
            [ch,Pop]=CheckRho(anyx,anyM,anyMs,anyNs,A);
        else
            
            [ch,Pop]=CheckRhoVlim(anyx,anyM,anyMs,anyNs);
        end
    end

thopt=0;
Popt=[];
Qopt=eye(n);
rhopt=0;


for it=1:length(th)
    V   =   vf*exp(1i*th(it));
    x   =   real(V);
    y   =   imag(V);
    M   =   sum(Matsx,3)*x+sum(Matsy,3)*y;
    
    
    for lt=1:10
%{
    [U,S,W] =   svd(M);
    Ms  =   Matsx(:,:,nsb);
    Ns  =   Matsy(:,:,nsb);
    M   =   S;  

    for jt=1:size(Ms,3)
        Ms(:,:,jt)  =   U'*Ms(:,:,jt)*W;
        Ns(:,:,jt)  =   U'*Ns(:,:,jt)*W;
    end
%}          
    rhos=   [rhopt,.5];
    for kt=1:5
        
        save(fname,'rhopt','Popt','Qopt');
        chk =   ffun(mean(rhos),M,Ms,Ns);
        if(chk)
            rhos=[mean(rhos),max(rhos)];
        else
            rhos=[min(rhos),mean(rhos)];
        end
    end
        
    rho=min(rhos);
    if(rho>rhopt)
        rhopt=rho;
        [~,P]=ffun(rhopt,M,Ms,Ns);
        Popt=P*U';
        Qopt=W;
        thopt=th(it);
    end
    disp(rhopt);
end



vref=vf*[cos(thopt);sin(thopt)];

end