function [Mv,Mbr,JF,nn]=MakePolarJac2(mpc)
    Y   =   makeYbus(mpc);
    n       =   size(mpc.bus,1);
sb      =   find(mpc.bus(:,2)==3);
nsb     =   setdiff(1:n,sb);
pv      =   find(mpc.bus(:,2)==2);
pq      =   find(mpc.bus(:,2)==1);
nsb     =   nsb(:);
pq      =   pq(:);
pv      =   pv(:);
npq     =   length(pq);
Vpv     =   mpc.bus(pv,8);Vpv=Vpv(:);
Vnom    =   mpc.bus(sb,8);Vnom=Vnom(:);
nbr     =   size(mpc.branch,1);

Y       =   makeYbus(mpc);
%Yd      =   exp(1i*(-pi/2+.01))./diag(Y);
%Y       =   diag(Yd)*Y;
Mv      =   cell(npq,1);
Mbr     =   cell(nbr,1);
for ii  =   1:npq
    [M,inds]=   Mati(pq(ii));    
    Mv{ii}  =   {M,inds,pq(ii)};
end

for ii  =   1:nbr
    ij          =   [mpc.branch(ii,1),mpc.branch(ii,2)];
    [Ms,Mc,inds]=   Matij(ij(1),ij(2));    
    Mbr{ii}  =   {Ms,Mc,inds,ij};
end


cc      =   randn(npq+n-1,1);
f       =   @(anyx) fg(anyx,cc);
xmin    =   minFunc(f,randn(npq+n-1,1),struct('derivativeCheck','on'));
JF      =   @Jacf;
nn      =   npq+n-1;

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

    function [s,c,v]=formscv(x)
        V   =   formV(x);
        s   =   imag(V(mpc.branch(:,1)).*conj(V(mpc.branch(:,2))));
        c   =   real(V(mpc.branch(:,1)).*conj(V(mpc.branch(:,2))));
        v   =   abs(V(pq)).^2;
    end

    function [f,grad]=fg(x,cc)
        f       =   cc'*getpf(x);     
        W       =   Jacf(x);
        grad    =   W'*cc;
    end

    function W=Jacf(x)
        [s,c,v] =   formscv(x);
        W       =   zeros(length(x));
        for it=1:npq
            Mvi =   Mv{it}{1};
            ind =   Mv{it}{2};
            W(ind,ind)  =   W(ind,ind)+Mvi*v(it);
        end
        
        for it=1:nbr
            Msi =   Mbr{it}{1};
            Mci =   Mbr{it}{2};
            ind =   Mbr{it}{3};            
            W(ind,ind)  =   W(ind,ind)+Msi*s(it)+Mci*c(it);
        end
    end

end