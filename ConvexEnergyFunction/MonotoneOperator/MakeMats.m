function [Matx,Maty]=MakeMats(mpc)
Y=makeYbus(mpc);
sb=find(mpc.bus(:,2)==3);
pv=find(mpc.bus(:,2)==2);
pq=find(mpc.bus(:,2)==1);
n=size(Y,1);
G=real(Y);
B=imag(Y);
Matsx=zeros(2*n,2*n,n);
Matsy=zeros(2*n,2*n,n);

for it=1:n
    e   =   sparse(it,1,1,n,1);
    Gi  =   G(it,:);
    Gip =   sparse(1,pq,G(it,pq),1,n);
    Bi  =   B(it,:);
    Bip =   sparse(1,pq,B(it,pq),1,n);
    if((mpc.bus(it,2)==1)||(mpc.bus(it,2)==3))    
        M   =   [diag(Gi)+e*Gi,diag(Bi)-e*Bi;...
                -diag(Bip)-e*Bi,diag(Gip)-e*Gi];
            
        N   =   [-diag(Bi)+e*Bi,diag(Gi)+e*Gi;...
                -diag(Gip)+e*Gi,-diag(Bip)-e*Bi];
    elseif(mpc.bus(it,2)==2)
        M   =   [diag(Gi)+e*Gi,diag(Bi)-e*Bi;...
                -diag(Bip)+2*e*e',diag(Gip)];
            
        N   =   [-diag(Bi)+e*Bi,diag(Gi)+e*Gi;...
                -diag(Gip),-diag(Bip)+2*e*e'];
    end
            
    Matsx(:,:,it)   =   M;
    Matsy(:,:,it)   =   N;
end



Matx=zeros(2*n-2,2*n-2,n);
Maty=zeros(2*n-2,2*n-2,n);

for it=1:n
    M   =   squeeze(Matsx(:,:,it));
    N   =   squeeze(Matsy(:,:,it));
    M([sb;n+sb],:) = [];
    M(:,[sb;n+sb]) = [];
    N([sb;n+sb],:) = [];
    N(:,[sb;n+sb]) = [];
    Matx(:,:,it)=M;
    Maty(:,:,it)=N;
end


end