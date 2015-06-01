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
%Yd      =   1./sqrt(abs(diag(Y)));
%Y       =   diag(Yd)*Y*diag(Yd);
Mv      =   cell(npq,1);
Mbr     =   cell(nbr,1);
for ii  =   1:npq
    [M,inds]=   Mati(pq(ii));    
    Mv{ii}  =   {M,inds,pq(ii)};
end

for ii  =   1:nbr
    [Ms,Mc,inds]=   Matij(mpc.branch(ii,1),mpc.branch(ii,2));    
    Mbr{ii}  =   {Ms,Mc,inds};
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

    function S=Sij(x,i,j)
        V   =   formV(x);
        S   =   V(i).*conj(Y(sub2ind(size(Y),i,j)).*V(j));
    end

    function f=pftwo(x)
        V   =   formV(x);
        S   =   V;
        for it=1:length(V)
            js          =   find(Y(it,:));
            js(js==it)  =   [];
            S(it)       =   sum(Sij(x,it(ones(length(js))),js));
        end
        f   =   [real(S(nsb));imag(S(pq))];
    end

    function [Ms,Mc,inds]=Matij(i,j)
        G   =   real(Y(i,j));
        B   =   imag(Y(i,j));
        Ms  =   [-G,G,B,B;...
                 -G,G,-B,-B;...
                 B,-B,G,G;...
                 B,-B,-G,-G];
        Mc  =   [B,-B,G,G;...
                 -B,B,G,G;...
                 G,-G,-B,-B;...
                 -G,G,-B,-B]; 
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