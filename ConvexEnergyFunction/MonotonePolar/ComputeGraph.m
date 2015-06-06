function Gn=ComputeGraph(mpc,Mv,Mbr,Qse,Qsi)
nbus    =   size(mpc.bus,1);
sb      =   find(mpc.bus(:,2)==3);
nsb     =   setdiff(1:nbus,sb);
pq      =   find(mpc.bus(:,2)==1);
nsb     =   nsb(:);
pq      =   pq(:);
npq     =   length(pq);
G       =   zeros(3*(nbus-1)+npq);
z       =   (1:(nbus-1+npq))';
V       =   zeros(nbus,2);
V(nsb,:)=   (npq+nbus-1)+[(1:((nbus-1)))',((nbus):2*(nbus-1))'];
for it=1:length(Mv)    
    ind     =   Mv{it}{2};
    ii      =   Mv{it}{3};
    indv    =   [z(ind);V(ii,:)'];
    indv(indv==0)   =   [];    
    G(indv,indv)    =   G(indv,indv)+1;
end

for it=1:length(Mbr)    
    ind     =   Mbr{it}{3};
    ij      =   Mbr{it}{4};
    indv    =   [z(ind);V(ij(1),:)';V(ij(2),:)'];
    indv(indv==0)   =   [];
    G(indv,indv)    =   G(indv,indv)+1;
end

for it=1:length(Qse)  
    indv    =   Qse{it}(:,1:2);    
    indv    =   unique(toind(indv(:)));
    
    indv(indv==0)   =   [];
    G(indv,indv)    =   G(indv,indv)+1;
end


for it=1:length(Qsi)  
    indv    =   Qsi{it}(:,1:2);
    indv    =   unique(toind(indv(:)));
    
    indv(indv==0)   =   [];
    G(indv,indv)    =   G(indv,indv)+1;
end

Gn   =   ComputeDecomposition(spones(G));

end
function res=toind(indv)
    res=zeros(length(indv),1);
    for it=1:length(res)
        res(it) =   toin(indv{it});
    end
end
function res=toin(anyx)
res=0;
if(iscell(anyx))
    if(ischar(anyx{1}))
        res=str2double(anyx);
    end
end
end