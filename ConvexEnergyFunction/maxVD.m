function delth=maxVD(mpc,delrho)
%define_constants;
Y       =   imag(makeYbus(mpc));
Yd      =   Y-diag(diag(Y));
n       =   size(Y,1);
m       =   size(mpc.branch,1);

PQ      =   (mpc.bus(:,2)==1);
m       =   (mpc.bus(:,2)==1)'*ones(n,1);

Yt      =   abs(Y(PQ,PQ));

Yss     =   Yt-diag(diag(Yt));

[i,j,ys]   =   find(Yss>0);
indi    =   find(i>j);
i       =   i(indi);
j       =   j(indi);
ys      =   ys(indi);
Ynn     =   abs(Yd(PQ,~PQ));

nk      =   length(i);
[ik,jk] =   find(Ynn>0);
mk      =   sum(sum(Ynn>0));

cvx_clear
cvx_begin
variable x(nk,1);
variable y(mk,1);
maximize(min([x;y]));
subject to
diag(2*diag(Yt)-delrho*sum(Ynn.*sparse(ik,jk,y,size(Ynn,1),size(Ynn,2)),2))-sparse(i,j,x.*ys,m,m)-sparse(i,j,x.*ys,m,m)'-sparse(i,i,ys.*x*delrho,m,m)-sparse(j,j,ys.*x*delrho,m,m)==semidefinite(m);
x>=1
y>=1
cvx_end
delth=acos(1./min([x;y]))*180/pi;
end 