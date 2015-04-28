function [Matsx,Matsw,Qssx,Qssw,fun,ff,TestQ,nn]=MakeMats(mpc,flim,sclV)
Y       =   makeYbus(mpc);
sb      =   find(mpc.bus(:,2)==3);
n       =   size(Y,1);
nsb     =   setdiff(1:n,sb);
G       =   real(Y);
B       =   imag(Y);
nbr     =   size(mpc.branch,1);
if(flim>0)
Amat    =   sparse(1:nbr,mpc.branch(:,1),1,nbr,n)...
                -sparse(1:nbr,mpc.branch(:,2),1,nbr,n);
else            
Amat    =   [];
nbr     =   0;
end

nn      =   4*(n-1)+nbr;
Vnom    =   exp(0*1i)*mpc.bus(sb,8);
Vmin    =   mpc.bus(nsb,13);
Vmax    =   mpc.bus(nsb,12);
Vavg    =   (Vmin+Vmax)/2;
Vmin    =   Vavg+sclV*(Vmin-Vavg);
Vmax    =   Vavg+sclV*(Vmax-Vavg);



            

Qss     =   MakeQs();
Mats    =   MakeMats();
Qssx    =   Qss(1:2*(n-1)+1,1:2*(n-1)+1,1:2*(n-1));
Qssw    =   Qss(1:2*(n-1)+1,1:2*(n-1)+1,2*n-1:end);
Matsx   =   Mats(1:2*(n-1),1:2*(n-1),1:2*n-1);
Matsw   =   Mats(1:2*(n-1),1:2*(n-1),2*n:end);

    function Mat=MakeMats()
        Mat         =   zeros(nn,nn,nn+1);
        Mat(:,:,1)  =   MakeJ(zeros(nn,1));
        for it=1:nn
            Mat(:,:,it+1)=MakeJ(sparse(it,1,1,nn,1))-Mat(:,:,1);
        end
    end
    function res=bsol(k)
        res=zeros(nn,1);
        if((k>=1)&&(k<=nn))
            res=sparse(k,1,1,nn,1);
        end
    end
    function Jac=MakeJ(V)
        Vx      =   real(Vnom)*ones(n,1);
        Vx(nsb) =   V(1:(n-1))+Vx(nsb);
        Vy      =   imag(Vnom)*ones(n,1);
        Vy(nsb) =   V(n:(2*n-2))+Vy(nsb);
        Gax     =   diag(G(nsb,:)*Vx);
        Gbx     =   diag(Vx(nsb))*G(nsb,nsb);
        Gay     =   diag(G(nsb,:)*Vy);
        Gby     =   diag(Vy(nsb))*G(nsb,nsb);
        Bax     =   diag(B(nsb,:)*Vx);
        Bbx     =   diag(Vx(nsb))*B(nsb,nsb);
        Bay     =   diag(B(nsb,:)*Vy);
        Bby     =   diag(Vy(nsb))*B(nsb,nsb);
        Jac     =   [Gax+Gbx+Bby-Bay,Gay+Gby+Bax-Bbx,zeros(n-1,2*(n-1)+nbr);...
                        Gby-Gay-Bbx-Bax,Gax-Gbx-Bby-Bay,zeros(n-1,2*(n-1)+nbr);...
                        2*diag(Vx(nsb)),2*diag(Vy(nsb)),-diag(V(2*n-1:3*n-3)),zeros((n-1),(n-1)+nbr);...
                        -2*diag(Vx(nsb)),-2*diag(Vy(nsb)),zeros(n-1,n-1),-diag(V(3*n-2:4*n-4)),zeros(n-1,nbr)];
        
        if(~isempty(Amat))
           Jac=[Jac;-2*diag(Amat*Vx)*Amat(:,nsb),-2*diag(Amat*Vy)*Amat(:,nsb),zeros(nbr,2*(n-1)),-diag(V(4*n-3:4*n-4+nbr));];
        end
    end

    function f=funf(V)
        Vx      =   real(Vnom)*ones(n,1);
        Vx(nsb) =   V(1:(n-1))+Vx(nsb);
        Vy      =   imag(Vnom)*ones(n,1);
        Vy(nsb) =   V(n:(2*n-2))+Vy(nsb);        
        Vc      =   Vx+1i*Vy;
        S       =   Vc.*conj(Y*Vc);
        f       =   [real(S(nsb));imag(S(nsb));...
            abs(Vc(nsb)).^2-V(2*n-1:3*n-3)-(Vmin).^2;...
             -abs(Vc(nsb)).^2-V(3*n-2:4*n-4)+(Vmax).^2];
         if(~isempty(Amat))
            f=[f;-abs(Amat*Vc).^2-V(4*n-3:4*n-4+nbr)+flim^2];
         end
    end

    function Qs=MakeQs
        Qs=zeros(nn+1,nn+1,nn);
        FF=funf(zeros(nn,1));
        for k=1:nn
            Qs(1,1,k)   =   FF(k);
        end
        for ii=2:(nn+1)
            FF=funf(bsol(ii-1));
            FFa=funf(2*bsol(ii-1));
            for k=1:nn
                Qs(ii,ii,k)   =   (FFa(k)-2*FF(k)+Qs(1,1,k))/2;
                Qs(1,ii,k)    =   (FF(k)-Qs(1,1,k)-Qs(ii,ii,k))/2;
                Qs(ii,1,k)    =   Qs(1,ii,k);
            end
        end
        for ii=2:(nn+1)
            for jj=(ii+1):(nn+1)
                FF=funf(bsol(ii-1)+bsol(jj-1));
                for k=1:nn
                    Qs(ii,jj,k)   =   (FF(k)-2*(Qs(1,ii,k)+Qs(1,jj,k))-Qs(1,1,k)-Qs(ii,ii,k)-Qs(jj,jj,k))/2;
                    Qs(jj,ii,k)   =   Qs(ii,jj,k);
                end
            end
        end
    end

    function [f,grad]=fres(V)
        V(2*n-1:end)=exp(V(2*n-1:end));
        F=funf(V);
        Jac=Mats(:,:,1)+sum(bsxfun(@times,Mats(:,:,2:end),permute(V,[3,2,1])),3);
        f=norm(F-1)^2;
        grad=2*Jac'*(F-1);
    end

    function chk=TestQs(V)
        Ft=zeros(nn,1);
        for kk=1:nn
            Ft(kk)=[1;V]'*Qss(:,:,kk)*[1;V];
        end
        chk=norm(Ft-funf(V));
    end
fun=@funf;
ff=@fres;
TestQ=@TestQs;
end