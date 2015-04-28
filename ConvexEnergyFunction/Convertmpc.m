efunction [Struct,mpc]=Convertmpc(mpc,alph,beta)

mpc.bus(:,4)   =    beta*mpc.bus(:,4);%.*(mpc.bus(:,4)>0)+(mpc.bus(:,4)<=0)*mean(mpc.bus(:,4)));
mpc.bus(:,3)   =    beta*mpc.bus(:,3);
mpc.gen(:,2)   =    beta*mpc.gen(:,2);


P           =   mpc.bus(:,3);
Q           =   mpc.bus(:,4);

P(mpc.gen(:,1))  =   P(mpc.gen(:,1))-mpc.gen(:,2);
if(abs(sum(P))>1e-5)
    return;
end
P           =   P/mpc.baseMVA;
Q           =   Q/mpc.baseMVA;
Qvar        =   alph*Q;
Q           =   Q-Qvar;
n           =   size(mpc.bus,1);
Bmat        =   sparse(mpc.branch(:,1),mpc.branch(:,2),1./mpc.branch(:,4),n,n);
Bmat        =   Bmat+Bmat';
PQInds      =   find(mpc.bus(:,2)==1);
PVInds      =   find((mpc.bus(:,2)==2)|(mpc.bus(:,2)==3));
Q           =   Q(PQInds);
Qvar        =   Qvar(PQInds);

c           =   min(Qvar./sum(Bmat(PQInds,PQInds),2));
Bmat        =   Bmat-diag(sum(Bmat,2));
Vs          =   mpc.bus(PVInds,8);

Struct=struct('P',P,'Q',Q,'Qvar',Qvar,'Bmat',Bmat,'PQInds',PQInds,'PVInds',PVInds,'Vs',Vs,'c',c);
end