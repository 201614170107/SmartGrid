function [chk,P]=CheckRhoVlim(rho,A,As,Bs)
Mxy=[1-sqrt(2),1;...
    -1,sqrt(2)-1;...
    -1,1-sqrt(2);...
    1-sqrt(2),-1];
Mxy=[Mxy;-Mxy(:,1),Mxy(:,2)];
%Mxy=[0,1;1,0;0,-1;-1,0];

Sym=@(anyx) (anyx+anyx');

n=size(A,1);
[ixx,jxx]=find(triu(true(n)));
[ix2,jx2]=find((true(n)));

SI=@(anyx)SymInds(ixx,jxx,anyx,n);
SI2=@(anyx) sparse(ix2,jx2,anyx,n,n);

k=size(As,3);
ni=length(ixx);
ni2=length(ix2);

cvx_clear;
cvx_precision('low');
cvx_begin quiet
variable P(n,n);
variable X(ni,k);
maximize(1)
subject to
for jt=1:k
    for kt=1:size(Mxy,1)
        SI(X(:,jt))-rho*(Mxy(kt,1)*Sym(P*As(:,:,jt)+Mxy(kt,2)*Sym(P*Bs(:,:,jt))))==semidefinite(n);
    end
end
Sym(P*A)-SI(sum(X,2))-1e-4*speye(n)==semidefinite(n)
cvx_end
chk=(cvx_optval>0.01);
end
