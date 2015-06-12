function [Mats,JF,nn]=MakePolarJac(mpc)
Y       =   makeYbus(mpc);
n       =   size(mpc.bus,1);
sb      =   find(mpc.bus(:,2)==3);
nsb     =   find(~(mpc.bus(:,2)==3));
pv      =   find(mpc.bus(:,2)==2);
pq      =   find(mpc.bus(:,2)==1);
nsb     =   nsb(:);
pq      =   pq(:);
pv      =   pv(:);
npq     =   length(pq);
Vpv     =   mpc.bus(pv,8);Vpv=Vpv(:);
Vnom    =   mpc.bus(sb,8);Vnom=Vnom(:);
nbr     =   size(mpc.branch,1);
bi      =   mpc.branch(:,1);
bj      =   mpc.branch(:,2);

Y       =   makeYbus(mpc);
vv      =   [1;arrayfun(@(anyx) {{anyx}},(n-1+npq)+(1:2*(n-1))')];
zz      =   (1:(n-1+npq))';
Mats    =   cell(npq+2*nbr,1);
for ii  =   1:npq
    [M,inds]    =   Mati(pq(ii)); 
    Q           =   GetQuads(@(anyx) Vi(anyx,pq(ii)),2*(n-1));
    Q           =   QuadPoly(Q,vv,3*(n-1)+npq);
    Mats{ii}    =   QuadMatPoly(M,zz(inds),Q);
end

for ii  =   1:nbr
    ij          =   [mpc.branch(ii,1),mpc.branch(ii,2)];
    [Ms,Mc,inds]=   Matij(ij(1),ij(2));    
    Qc          =   GetQuads(@(anyx) Vijc(anyx,ij(1),ij(2)),2*(n-1));
    Q           =   QuadPoly(Qc,vv,3*(n-1)+npq);
    Mats{npq+2*(ii-1)+1}  =   QuadMatPoly(Ms,zz(inds),Q);
    Qc          =   GetQuads(@(anyx) Vijr(anyx,ij(1),ij(2)),2*(n-1));
    Q           =   QuadPoly(Qc,vv,3*(n-1)+npq);
    Mats{npq+2*(ii-1)+2}  =   QuadMatPoly(Mc,zz(inds),Q);
end
x       =   randn(2*(n-1),1);
Vx    =   formVc(x);

cc      =   randn(npq+n-1,1);
f       =   @(anyx) fg(anyx,cc);
xmin    =   minFunc(f,randn(npq+n-1,1),struct('derivativeCheck','on'));
JF      =   @Jacf;
nn      =   npq+n-1;

    function V=formVc(x)
        V       =   zeros(n,1);
        V(sb)   =   Vnom;
        V(nsb)  =   x(1:(n-1))+1i*x(n:end);
    end
    
     function res=Vi(x,i)
        V       =   formVc(x);
        res     =   abs(V(i)).^2;
     end

    function res=Vijr(x,i,j)
        V       =   formVc(x);
        res     =   real(V(i).*conj(V(j)));
    end

    function res=Vijc(x,i,j)
        V       =   formVc(x);
        res     =   imag(V(i).*conj(V(j)));
    end

    function J=MakeJac(V)
        J   =   zeros(n+npq-1);
        for it=1:npq
             [Mit,ins]      =   Mati(pq(it)); 
             J(ins,ins)     =   J(ins,ins)+Mit*abs(V(pq(it))).^2;
        end
        for it=1:nbr
             ijc            =   [bi(it),bj(it)];
             [Ma,Mb,ins]    =   Matij(ijc(1),ijc(2)); 
             J(ins,ins)     =   J(ins,ins)+Ma*imag(V(ijc(1))*conj(V(ijc(2))));
             J(ins,ins)     =   J(ins,ins)+Mb*real(V(ijc(1))*conj(V(ijc(2))));
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

   
    function [Ms,Mc,inds]=Matij(i,j)
        Gi  =   real(Y(i,j));
        Bi  =   imag(Y(i,j));
        Gj  =   real(Y(j,i));
        Bj  =   imag(Y(j,i));
        Ms  =   [-Gi,Gi,Bi,Bi;...
                 -Gj,Gj,-Bj,-Bj;...
                 Bi,-Bi,Gi,Gi;...
                 Bj,-Bj,-Gj,-Gj];
        Mc  =   [Bi,-Bi,Gi,Gi;...
                 -Bj,Bj,Gj,Gj;...
                 Gi,-Gi,-Bi,-Bi;...
                 -Gj,Gj,-Bj,-Bj]; 
        [Ms,~]      =   shrinkM(Ms,[i;j;n+i;n+j]);     
        [Mc,inds]   =   shrinkM(Mc,[i;j;n+i;n+j]);    
    end

    function [M,inds]=shrinkM(Ms,inds)        
        M   =   zeros(2*n);
        M(inds,inds)    =   Ms;
        M   =   M([nsb;n+pq],[nsb;n+pq]);
        inds=   find(sum(abs(M),2));
        M   =   M(inds,inds);
    end

    function [Mx,inds]=Mati(i)
        G   =   real(Y(i,i));
        B   =   imag(Y(i,i));
        Mx  =   2*[0,G;0,-B];
        [Mx,inds]   =   shrinkM(Mx,[i;n+i]);
    end

    function [f,grad]=fg(x,cc)
        f       =   cc'*getpf(x);     
        W       =   MakeJac(formV(x));
        grad    =   W'*cc;
    end
end