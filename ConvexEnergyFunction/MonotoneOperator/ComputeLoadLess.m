function VSr=ComputeLoadLess(mpc)
define_constants;
sb=find((mpc.bus(:,2)==3));
Y=makeYbus(mpc);
nb=size(Y,1);
nsb=setdiff(1:nb,sb);
Y=Y(nsb,:);
VS=zeros(length(nsb));
Vsb=mpc.bus(sb,VM)*exp(1i*mpc.bus(sb,VA)*pi/180);


nbrn=size(mpc.branch,1);
Amat=sparse(mpc.branch(:,1),1:nbrn,1,nb,nbrn)-sparse(mpc.branch(:,2),1:nbrn,1,nb,nbrn);
nsb=find(~(mpc.bus(:,2)==3));
A   =   Amat(nsb,:);

for i=1:length(nsb)
    Mi      =   Y;
    Mi(i,:) =   sparse(1,nsb(i),1,1,nb);
    bi      =   -Y(:,sb)*Vsb;
    bi(i)   =   0;
    VS(:,i) =   (Mi(:,nsb))\bi;
end

VSr =   zeros(size(mpc.bus,1),length(nsb));
VSr(nsb,:)   =   VS;
VSr(sb,:)    =   Vsb;
end