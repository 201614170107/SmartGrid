function [Mats,Qse,Qsi,nn,ff,nb,JF,bds,indpairs]=MakePolarJac(mpc,del)
n       =   size(mpc.bus,1);
sb      =   find(mpc.bus(:,2)==3);
nsb     =   setdiff(1:n,sb);
pv      =   find(mpc.bus(:,2)==2);
pq      =   setdiff(nsb,pv);
nsb     =   nsb(:);
pq      =   pq(:);
pv      =   pv(:);
npq     =   length(pq);
Vpv     =   mpc.bus(pv,8);Vpv=Vpv(:);
Vnom    =   mpc.bus(sb,8);Vnom=Vnom(:);
Vmin    =   mpc.bus(:,13);Vmin=Vmin(:);
Vmax    =   mpc.bus(:,12);Vmax=Vmax(:);
Vmin(pv)=   Vpv;
Vmax(pv)=   Vpv;
nbr     =   size(mpc.branch,1);
A       =   sparse(1:nbr,mpc.branch(:,1),1,nbr,n)...
                -sparse(1:nbr,mpc.branch(:,2),1,nbr,n);
A       =   A';
ijs     =   mpc.branch(:,1:2);

Y       =   makeYbus(mpc);
Gd      =   real(diag(Y));
Bd      =   imag(diag(Y));
Ybr     =   Y(sub2ind(size(Y),mpc.branch(:,1),mpc.branch(:,2)));
G       =   real(Ybr);
B       =   imag(Ybr);
cc      =   randn(length(nsb)+length(pq),1);
nn      =   length(cc);
[~,fi]  =   func(zeros(2*nbr+npq));
Qs      =   GetQuads(@func,2*nbr+npq);
Qse     =   Qs(:,:,fi==1);
Qsi     =   Qs(:,:,fi==0);
Mats    =   MakeMats();
JF      =   @Jac;
nb      =   2*nbr+npq;
bds     =   [Vmin(pq).^2,Vmax(pq).^2]-1;
bds     =   [bds;[Vmin(ijs(:,1)).*Vmin(ijs(:,2))*cos(del)-1,Vmax(ijs(:,1)).*Vmax(ijs(:,2))-1]];
bds     =   [bds;[-Vmax(ijs(:,1)).*Vmax(ijs(:,2))*sin(del),Vmax(ijs(:,1)).*Vmax(ijs(:,2))*sin(del)]];
indpairs=   [mpc.branch(:,1),mpc.branch(:,2),((npq+1):(npq+nbr))',((npq+nbr+1):(npq+2*nbr))'];
%Matpairs=   cell(nbr,


    function Mats=MakeMats()
        Mats=zeros(length(nsb)+npq,length(nsb)+npq,2*nbr+npq+1);
        Mats(:,:,1) =   Jac(zeros(2*nbr+npq,1));
        for it=1:2*nbr+npq
            Mats(:,:,it+1)  =   Jac(sparse(it,1,1,2*nbr+npq,1))-...
                                    Mats(:,:,1);
        end
    end
    function V=formV(x)
        V       =   zeros(n,1);
        V(sb)   =   Vnom;
        V(nsb)  =   exp(1i*x(1:n-1));
        V(pv)   =   V(pv).*Vpv;
        V(pq)   =   V(pq).*exp(x(n:end));
    end

    function f=getpf(x)
        V   =   formV(x);
        S   =   V.*conj(Y*V);
        f   =   [real(S(nsb));imag(S(pq))];
    end

    function J=Jac(x)
        Vpq     =   x(1:npq)+1;
        c       =   x((npq+1):(npq+nbr))+1;
        s       =   x((npq+nbr+1):(npq+2*nbr));
        Vs      =   zeros(n,1);
        Vs(pq)  =   Vpq;
        Vs(pv)  =   Vpv.^2;
        Vs(sb)  =   Vnom.^2;
        J       =   [A*diag(B.*c)*A'-abs(A)*diag(G.*s)*A',...
            2*diag(Gd.*Vs)+abs(A)*diag(G.*c)*abs(A)'+A*diag(B.*s)*abs(A)'];
        J       =   [J;...
                        A*diag(G.*c)*A'+abs(A)*diag(B.*s)*A',...
                        -2*diag(Bd.*Vs)+A*diag(G.*s)*abs(A)'+...
                        -abs(A)*diag(B.*c)*abs(A)'];
        J       =   J([nsb;n+pq],[nsb;n+pq]);
    end

    function [f,fi]=func(x)
        Vpq     =   x(1:npq);
        Vpq     =   Vpq(:);
        V       =   sparse([sb;pv;pq],1,[Vnom.^2;Vpv.^2;Vpq],n,1);
        c       =   x(((npq+1):(npq+nbr))');
        s       =   x(((npq+nbr+1):(npq+2*nbr))');
        f       =   c.^2+s.^2-V(ijs(:,1)).*V(ijs(:,2));      
        g       =   [c-del*V(ijs(:,1));...
                        c-del*V(ijs(:,2));Vpq-Vmin(pq).^2;Vmax(pq).^2-Vpq];
        g       =   [g;(c-del*V(ijs(:,1))).*(V(ijs(:,1))-Vmin(ijs(:,1)).^2)];
        g       =   [g;(c-del*V(ijs(:,1))).*(V(ijs(:,2))-Vmin(ijs(:,2)).^2)];
        g       =   [g;(c-del*V(ijs(:,2))).*(V(ijs(:,2))-Vmin(ijs(:,2)).^2)];
        g       =   [g;(c-del*V(ijs(:,2))).*(V(ijs(:,1))-Vmin(ijs(:,1)).^2)];

      %  g       =   g*g';
      %  g       =   g(triu(true(size(g,1))));        
        fi      =   [ones(size(f));zeros(size(g))];
        f       =   [f;g];
    end

    function [f,grad]=formf(x)
        V       =   sparse([sb;pv;pq],1,[Vnom;Vpv;exp(x(n:end))],n,1);
        th      =   sparse(nsb,1,x(1:n-1),n,1);
        c       =   exp(abs(A)'*log(V)).*cos(A'*th);
        s       =   exp(abs(A)'*log(V)).*sin(A'*th);
        J       =   Jac([exp(2*x(n:end))-1;c-1;s]);
        f       =   getpf(x);
        f2      =   abs(A)*diag(G)*c+A*diag(B)*s+Gd.*V.^2;
        f2      =   [f2;-abs(A)*diag(B)*c+A*diag(G)*s-Bd.*V.^2];
        f2      =   [f2(nsb);f2(n+pq)];
        f       =   cc'*f;        
        grad    =   J'*cc;
    end

ff=@formf;
end