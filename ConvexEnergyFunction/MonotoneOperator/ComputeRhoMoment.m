function [chk,sol]=ComputeRhoMoment(Mx,Qse,Qsw,di,usecvx)
nn          =   size(Mx,1);

if(di)
    Mscl        =   eye(nn);diag(diag(Mx(:,:,1)\eye(nn))); 
else
    Mscl        =   (Mx(:,:,1)\eye(nn)); 
end
for it=1:size(Mx,3)
    Mx(:,:,it)=Mscl*Mx(:,:,it);
end

Qsee        =   cell(size(Qse,3),1);
Qsew        =   cell(size(Qsw,3),1);
Qsobj       =   cell(nn,1);

sp          =   @(anyn) sparse([],[],[],anyn,anyn);

for it=1:size(Qse,3)
    Qsee{it}        =   blkdiag(squeeze(Qse(:,:,it)),sp(nn));    
end

for it=1:nn
    Mc              =   squeeze(Mx(it,:,:));
    Qsobj{it}       =   [sp(nn+1),Mc';Mc,sp(nn)];    
end


for it=1:size(Qsw,3)
    Qsew{it}        =   blkdiag(squeeze(Qsw(:,:,it)),sp(nn));   
end



Qsee        =   Cleanup(Qsee);
Qsew        =   Cleanup(Qsew);
Qsobj       =   Cleanup(Qsobj);
QQ          =   sparse([1,nn+2:(2*nn+1)],[1,nn+2:(2*nn+1)],[1;-ones(nn,1)],2*nn+1,2*nn+1);

[chk,sol]   =   MomentRelaxation(Qsee,Qsew,Qsobj,QQ,usecvx);
end

function Qs=Cleanup(Qs)
    for it=1:length(Qs)
        Qs{it}(abs(Qs{it})<1e-5*max(abs(Qs{it}(:))))=0;
        Qs{it}=sparse(Qs{it});
    end
end
