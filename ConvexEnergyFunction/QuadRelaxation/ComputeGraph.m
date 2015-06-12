function Gn=ComputeGraph(Mats,Qs)
n       =   Qs{1}.n;
G       =   zeros(n);
for it=1:length(Mats)
    ind             =   Mats{it}.zi;
    is              =   Mats{it}.Qp.is;
    js              =   Mats{it}.Qp.js;
    for jt=1:length(is)
        i               =   toind(is{jt});
        j               =   toind(js{jt});
        indv            =   unique([ind;i;j]);
        indv(indv==0)   =   [];
        G(indv,indv)    =   G(indv,indv)+1;        
    end
end

for it=1:length(Qs)
    indv            =   Qs{it}.GetInds();
    G(indv,indv)    =   G(indv,indv)+1;
end

Gn   =   ComputeDecomposition(spones(G));

end
function res=toind(i)
res =   0;
if(ischar(i))
    res =   str2double(i);
end
end