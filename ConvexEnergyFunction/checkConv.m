function p=checkConv(rho,theta,mpc)
define_constants;
Y       =   imag(makeYbus(mpc));
n       =   size(Y,1);
m       =   size(mpc.branch,1);

PQ      =   (mpc.bus(:,BUS_TYPE)==1);
PV      =   ~(PQ);
A       =    sparse(1:m,mpc.branch(:,1),1,m,n)...
                -sparse(1:m,mpc.branch(:,2),1,m,n);
if(any(abs(A*theta)>pi/2))
    p=0;    
end

Ymat    =   sparse(mpc.branch(:,1),mpc.branch(:,2),...
                        (2-exp(A*rho)./cos(A*theta))./mpc.branch(:,4),n,n);
Ymat    =   (Ymat+Ymat');                    
Y       =   diag(sum(Ymat,2));
Ymat    =   sparse(mpc.branch(:,1),mpc.branch(:,2),...
                        (1./cos(A*theta))./mpc.branch(:,4),n,n);
Y       =   Y-Ymat-Ymat';

[~,p]   =   chol(Y(PQ,PQ));
p       =   (p==0);
end