function chk=CheckEF(mpc,Vs)
    Y=makeYbus(mpc);
    B=-imag(Y);
    n=size(B,1);
    pq=(mpc.bus(:,2)==1);
    VM=abs(Vs);
    VA=Vs./VM;
    i=mpc.branch(:,1);
    j=mpc.branch(:,2);
    M=sparse(i,j,(VM(i)./VM(j)).^2./(real(VA(i)./VA(j))),n,n).*B;
    M=M+M';
    M=M-diag(sum(M,2));
    [~,chk]=chol(M(pq,pq));
    chk=(chk==0);
end