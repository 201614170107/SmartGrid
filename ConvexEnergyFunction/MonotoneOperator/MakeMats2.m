function [Mats,Qse,Qsi,nn,TQs,ff]=MakeMats2(mpc,flim,sclV)
Y       =   makeYbus(mpc);
sb      =   find(mpc.bus(:,2)==3);
n       =   size(Y,1);
nsb     =   setdiff(1:n,sb);
pv      =   find(mpc.bus(:,2)==2);
pq      =   setdiff(nsb,pv);
G       =   real(Y);
B       =   imag(Y);
nbr     =   size(mpc.branch,1);
Ai      =   [];
Aj      =   [];
if(flim>0)
    Amat    =   sparse(1:nbr,mpc.branch(:,1),1,nbr,n)...
        -sparse(1:nbr,mpc.branch(:,2),1,nbr,n);
    Ai      =   mpc.branch(:,1);
    Aj      =   mpc.branch(:,2);
else
    Amat    =   [];
    nbr     =   0;
end

nn      =   2*(n-1);
Vnom    =   exp(0*1i)*mpc.bus(sb,8);
Vmin    =   mpc.bus(:,13);
Vmax    =   mpc.bus(:,12);
Vavg    =   (Vmin+Vmax)/2;
Vmin    =   Vavg+sclV*(Vmin-Vavg);
Vmax    =   Vavg+sclV*(Vmax-Vavg);
Vpv     =   mpc.bus(pv,8);


[~,Jac,Mats]    =   GetQuads(@funf,nn);
[Qse,Qsi]       =   MakeQs();
Qse(abs(Qse)<1e-10)  =   0;
Qsi(abs(Qsi)<1e-10)  =   0;


    %   Power Flow Operator %
    function f=funf(V)
        [Vx,Vy] =   formV(V);
        Vc      =   Vx+1i*Vy;
        S       =   Vc.*conj(Y*Vc);
        f       =   [real(S);imag(S)];
        f(n+pv) =   abs(Vc(pv)).^2-Vpv.^2;
        f       =   [f(nsb);f(n+nsb)];
    end

    %   Helper Function - Add Slack Bus Voltage Reference %
    function [Vx,Vy,A]=formV(V)
        
        A       =   eye(2*n);
        A       =   A(:,[nsb;n+nsb]);
        Vx      =   A(1:n,:)*V+real(Vnom);
        Vy      =   A(n+1:2*n,:)*V+imag(Vnom);
    end

    %   Constraints on Voltages %
    function f=func(V)
        [Vx,Vy] =   formV(V);
        Vc      =   Vx+1i*Vy;
        f       =   [abs(Vc(pv)).^2-Vpv.^2;abs(Vc(pq)).^2-Vmin(pq).^2;...
            Vmax(pq).^2-abs(Vc(pq)).^2];
        if(~isempty(Amat))            
           % f=[f;real(Vc(Ai).*conj(Vc(Aj)))-flim*abs(Vc(Ai)).^2;real(Vc(Ai).*conj(Vc(Aj)))-flim*abs(Vc(Aj)).^2];
            f=[f;-abs(Amat*Vc).^2+flim^2];
        end
    end

    %   Constraints on Voltages in Quadratic Form   %
    function [Qseq,Qsiq]=MakeQs
        Qs      =   GetQuads(@func,nn);
        Qseq    =   Qs(:,:,1:length(pv));
        Qsiq    =   Qs(:,:,length(pv)+1:end);
    end

    %   Function to test correctness of Jacobian   %
    function [f,grad]=fres(V)
        F=funf(V);
        Jac=Mats(:,:,1)+sum(bsxfun(@times,Mats(:,:,2:end),permute(V,[3,2,1])),3);
        f=sum(F);
        grad=sum(Jac,1)';
    end

    %   Function to test correctness of Quadratic Form Constraints   %
    function chk=TestQs(V)
        Ft=func(zeros(nn,1));
        Qss=cat(3,Qse,Qsi);
        for kk=1:length(Ft)
            Ft(kk)=[1;V]'*Qss(:,:,kk)*[1;V];
        end
        chk=norm(Ft-func(V));
    end
ff=@fres;
TQs=@TestQs;
end