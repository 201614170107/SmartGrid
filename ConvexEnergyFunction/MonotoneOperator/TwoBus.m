clearvars;
clc;
mpc=loadcase('case2');
Y=makeYbus(mpc);

[Matx,Maty]=MakeMats(mpc);
Ma=sum(Matx,3);
Mb=sum(Maty,3);
M=Matx(:,:,2);
N=Maty(:,:,2);

Sym=@(anyx) (anyx+anyx');

xo=1;
yo=0;
Mo=Ma*xo+Mb*yo;

sols=[0,0];
sol2=-Y(1,2)*(xo+1i*yo)/Y(2,2);
sols=[sols;real(sol2),imag(sol2)];
sols=sols-[xo,yo;xo,yo];

%%

Xa=[-1,-.51;1,-.51;-1,-1.5;1,-1.5];
Xa=[Xa(:,2),Xa(:,1)];
cvx_begin
variable Wa(2,2)
minimize(1)
subject to
for it=1:size(Xa,1)
    Sym(Wa*(Mo+M*Xa(it,1)+N*Xa(it,2)))-1e-3*eye(2)==semidefinite(2);
end
cvx_end

%%

Xb=[-1,-.49;1,-.49;-1,.5;1,.5];
Xb=[Xb(:,2),Xb(:,1)];
cvx_begin
variable Wb(2,2)
minimize(1)
subject to
for it=1:size(Xb,1)
    Sym(Wb*(Mo+M*Xb(it,1)+N*Xb(it,2)))-1e-3*eye(2)==semidefinite(2);
end
cvx_end




%%
[xs,ys]=ndgrid(linspace(-1.5,.5,100),linspace(-1,1,100));
X=[xs(:),ys(:)];
figure(1);clf;

W=Wa;
chk=zeros(size(X,1),1);
for it=1:size(X,1)
    [~,p]=chol(Sym(W*(Mo+M*X(it,1)+N*X(it,2))));
    chk(it)=(p==0);
end

Y=X(chk>0,:);

figure(1);hold on;
scatter(Y(:,1),Y(:,2),20,'c');
xlim([-1 1]);
ylim([-1.5 .5]);


W=Wb;
chk=zeros(size(X,1),1);
for it=1:size(X,1)
    [~,p]=chol(Sym(W*(Mo+M*X(it,1)+N*X(it,2))));
    chk(it)=(p==0);
end

Y=X(chk>0,:);

figure(1);hold on;
scatter(Y(:,1),Y(:,2),20,'g');
ylim([-1 1]);
xlim([-1.5 .5]);

scatter(sols(:,1),sols(:,2),100,'r','filled');

