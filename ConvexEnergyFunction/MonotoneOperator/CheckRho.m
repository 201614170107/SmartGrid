function [chk,P]=CheckRho(rho,A,As,Bs,Z)
%Mxy=[1-sqrt(2),1;...
 %   -1,sqrt(2)-1;...
  %  -1,1-sqrt(2);...
   % 1-sqrt(2),-1];
%Mxy=[Mxy;-Mxy(:,1),Mxy(:,2)];
Mxy=[0,1;1,0;0,-1;-1,0];

    
Sym=@(anyx) (anyx+anyx');
n=size(A,1);
[ixx,jxx]=find(triu(ones(n)));
[ix2,jx2]=find((ones(n)));
ni=length(ixx);
ni2=length(ix2);

SI=@(anyx)SymInds(ixx,jxx,anyx,n);
SI2=@(anyx) anyx;%sparse(ix2,jx2,anyx,n,n);


[m,k]=size(Z);
if(~((m==size(As,3))&&(m==size(Bs,3))))
     ME = MException('Sizes of Inputs Incorrect');
    throw(ME);
end

cvx_clear;
cvx_quiet;
cvx_precision('low');
cvx_begin
variable P(n,n);
variable X(ni,k);
variable T(ni,k);
variable S(ni,k);
maximize(1)
subject to
for jt=1:k
    for kt=1:size(Mxy,1)
        SI(X(:,jt))-rho*(Mxy(kt,1)*SI(T(:,jt))+Mxy(kt,2)*SI(S(:,jt)))==semidefinite(n);              
    end
end
for jt=1:m
    Sym(SI2(P)*As(:,:,jt))==SI(T*Z(jt,:)');
    Sym(SI2(P)*Bs(:,:,jt))==SI(S*Z(jt,:)');
end
Sym(SI2(P)*A)-SI(sum(X,2))-1e-4*speye(n)==semidefinite(n)
cvx_end
chk=(cvx_optval>0.01);
P=SI2(P);
end