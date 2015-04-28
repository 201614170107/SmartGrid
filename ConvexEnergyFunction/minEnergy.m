function minfE=minEnergy(mpc)
    sb  =   find(mpc.bus(:,2)==3);
    n   =   size(mpc.bus,1);
    m   =   size(mpc.branch,1);
    
    PV  =   (mpc.bus(:,2)==2)|(mpc.bus(:,2)==1);
    PQ  =   mpc.bus(:,2)==1;
    P   =   -mpc.bus(:,3);
    Q   =   -mpc.bus(:,4);
    P(mpc.gen(:,1)) =   P(mpc.gen(:,1))+mpc.gen(:,2);
    Q(mpc.gen(:,1)) =   Q(mpc.gen(:,1))+mpc.gen(:,3);
    P   =   P/mpc.baseMVA;
    Q   =   Q/mpc.baseMVA;
    
    
    A   =   sparse(1:m,mpc.branch(:,1),1,m,n)...
                -sparse(1:m,mpc.branch(:,2),1,m,n);
    Ap  =   sparse(1:m,mpc.branch(:,1),1,m,n)...
                +sparse(1:m,mpc.branch(:,2),1,m,n);
    Bs  =   1./mpc.branch(:,4);
    
    Bm  =   sparse(mpc.branch(:,1),mpc.branch(:,2),Bs,n,n);
    Bm  =   Bm+Bm';
    Bi  =   sum(Bm,2);
    
    
    function [f,grad]=minE(rhoth)
        rho         =   rhoth(1:sum(PQ));
        theta       =   rhoth(sum(PQ)+1:end);
        th          =   zeros(n,1);
        rh          =   zeros(n,1);
        th(PV)      =   theta;
        rh(PQ)      =   rho;
        rh(~PQ)     =   log(mpc.bus(~PQ,8));
        
        p           =   checkConv(rh,th,mpc);
        if(~p)
            f   =   Inf;
            grad=   ones(size(rhoth));
            return;
        end
        
        
        f           =   -P'*th-Q'*rh-sum(exp(Ap*rh).*cos(A*th).*Bs)+.5*Bi'*exp(2*rh);
        gth         =   -P+A'*(exp(Ap*rh).*sin(A*th).*Bs);
        grh         =   -Q-Ap'*(exp(Ap*rh).*cos(A*th).*Bs)+Bi.*exp(2*rh);
        
        grad        =   [grh(PQ);gth(PV)];
    end

minfE=@minE;
end