function [chk,P]=CheckRhoLP(rho,A,As,Z)
Mxy=[1-sqrt(2),1;...
    -1,sqrt(2)-1;...
    -1,1-sqrt(2);...
    1-sqrt(2),-1];
Mxy=[Mxy;-Mxy(:,1),Mxy(:,2)];
%Mxy=[1,0;-1,0;0,1;0,-1];
    

n=size(As,1);
[ixx,jxx]=find(triu(abs(A)));
ni=length(ixx);
k=size(Z,2);
SI=@(anyx)SymInds(anyx,ixx,jxx,n);
[ix2,jx2]=find(triu(abs(A)));
SI2=@(anyx) sparse(ix2,jx2,anyx,n,n);
ni2=length(ix2);


Z=full(Z);
Z=full(Z);
m=size(Z,1);

cvx_clear;
cvx_quiet;
cvx_begin
variable P(ni2,1);
variable X(ni,k);
variable T(ni,k);
variable S(ni,k);
variable M(ni,k,size(Mxy,1));
variable N(ni,k,size(Mxy,1));
variable O(ni,k,size(Mxy,1));
variable Ml(ni,1);
variable Nl(ni,1);
variable Ol(ni,1);
maximize(1)
subject to
for jt=1:k
    for kt=1:size(Mxy,1)
        SI(X(:,jt))-rho*(Mxy(kt,1)*SI(T(:,jt))+Mxy(kt,2)*SI(S(:,jt)))...
            ==sparse(ixx,ixx,M(:,jt,kt),n,n)...
                +sparse(jxx,jxx,N(:,jt,kt),n,n)...
                +(sparse(ixx,jxx,O(:,jt,kt),n,n)+sparse(ixx,jxx,O(:,jt,kt),n,n)');
        for lt=1:ni
            norm([M(lt,jt,kt)-N(lt,jt,kt);2*O(lt,jt,kt)])<=(M(lt,jt,kt)+N(lt,jt,kt));
        end
    end
    
end
for jt=1:m
    As(:,:,jt)'*SI2(P)'+SI2(P)*As(:,:,jt)==SI(T*Z(jt,:)');
    As(:,:,m+jt)'*SI2(P)'+SI2(P)*As(:,:,m+jt)==SI(S*Z(jt,:)');
end
A'*SI2(P)'+SI2(P)*A-SI(sum(X,2))-1e-2*speye(n)==sparse(ix2,ix2,Ml,n,n)...
                +sparse(jx2,jx2,Nl,n,n)...
                +(sparse(ix2,jx2,Ol,n,n)+sparse(ix2,jx2,Ol,n,n)');
for lt=1:ni
    norm([Ml(lt)-Nl(lt);2*Ol(lt)])<=(Ml(lt)+Nl(lt));
end
cvx_end
chk=(cvx_optval>0.01);
P=SI2(P);
end