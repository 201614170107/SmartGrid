function rhopt=MaxRho(mpc,P,vref)
[As,Bs]=MakeMats(mpc);
sb=find(mpc.bus(:,2)==3);
nsb=setdiff(1:size(mpc.bus,1),sb);
tp  =   mpc.bus(:,2);
vs  =   mpc.bus(:,8);
tp  =   tp(nsb);
vs  =   -vs(nsb).^2+norm(vref)^2;
pq  =   find(tp==1);
pv  =   find(tp==2);

%Z   =   sparse(mpc.branch(:,1),1:size(mpc.branch,1),1,size(mpc.bus,1),size(mpc.branch,1))-...
 %           -sparse(mpc.branch(:,2),1:size(mpc.branch,1),1,size(mpc.bus,1),size(mpc.branch,1));
%Z   =   Z(nsb,:);        



A   =   sum(As,3)*vref(1)+sum(Bs,3)*vref(2);
As  =   As(:,:,nsb);
Bs  =   Bs(:,:,nsb);

Mxy=[1-sqrt(2),1;...
    -1,sqrt(2)-1;...
    -1,1-sqrt(2);...
    1-sqrt(2),-1];
Mxy=[Mxy;-Mxy(:,1),Mxy(:,2)];
%Mxy=[0,1;1,0;0,-1;-1,0];


n=size(A,1);
Sym=@(anyx) (anyx+anyx');

A=Sym(P*A);
for it=1:size(As,3)
    As(:,:,it)=Sym(P*As(:,:,it));
    Bs(:,:,it)=Sym(P*Bs(:,:,it));
end
k=size(As,3);

[ix,jx]=find(triu(ones(n)));%eye(n)+sum(abs(As),3)+sum(abs(Bs),3)));
ni=length(ix);
SI=@(anyx) SymInds(ix,jx,anyx,n);

rhos=[0,.5];
rhopt=0;
for itout=1:10
    rho=mean(rhos);
    cvx_begin
    variable X(ni,k);
    variable T(ni,length(pv));
    maximize(1)
    subject to
    for jtt=1:length(pq)
        jt=pq(jtt);
        for kt=1:size(Mxy,1)
            SI(X(:,jt))-rho*(Mxy(kt,1)*As(:,:,jt)+Mxy(kt,2)*Bs(:,:,jt))==semidefinite(n);
        end
    end
    for jtt=1:length(pv)
        jt=pv(jtt);
        SI(T(:,jtt))==semidefinite(n);
        for kt=1:size(Mxy,1)
            SI(X(:,jt))+SI(T(:,jtt))*(rho+vs(jt))--rho*(+Mxy(kt,1)*(As(:,:,jt)-SI(T(:,jtt)))+Mxy(kt,2)*(Bs(:,:,jt)-SI(T(:,jtt))))==semidefinite(n);
        end
    end
    A-SI(sum(X,2))==semidefinite(n)
    cvx_end
    chk=(cvx_optval>0.01);
    if(chk)
        rhopt=rho;
        rhos=[rho,rhos(2)];
    else
        rhos=[rhos(1),rho];
    end
end
end