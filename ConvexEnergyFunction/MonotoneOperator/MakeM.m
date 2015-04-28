function M=MakeM(mpc,x,y)
Y=-makeYbus(mpc);
sb=find(mpc.bus(:,2)==3);
pv=find(mpc.bus(:,2)==2);
n=size(Y,1);
G=real(Y);
Gd=-diag(G);
G(eye(3)>0)=0;
B=-imag(Y);
Bd=-diag(B);
B(eye(3)>0)=0;

    Ma              =   -diag(x)*G-diag(y)*B;
    Ma(eye(n)>0)    =   2*Gd.*x-G*x+B*y;
    Mb              =   -diag(y)*G+diag(x)*B;
    Mb(eye(n)>0)    =   2*Gd.*y-G*y-B*x;
    Mc              =   -diag(x)*B+diag(y)*G;
    Mc(eye(n)>0)    =   2*Bd.*x-B*x-G*y;
    Md              =   -diag(y)*B-diag(x)*G;
    Md(eye(n)>0)    =   2*Bd.*y-B*y+G*x;    
    Mc(pv,pv)       =   2*diag(x(pv));
    Md(pv,pv)       =   2*diag(y(pv));
    
    ns              =   setdiff(1:n,sb);
    M               =   [Ma(ns,ns),Mb(ns,ns);Mc(ns,ns),Md(ns,ns)];
end
