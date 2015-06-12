function [Qse,Qsi,Qss,nn,TQs,ff]=MakeMats(mpc,gam,plims,qlims)
Y       =   makeYbus(mpc);
sb      =   find(mpc.bus(:,2)==3);
n       =   size(Y,1);
nsb     =   setdiff(1:n,sb);
pv      =   find(mpc.bus(:,2)==2);
pq      =   setdiff(nsb,pv);
nn      =   2*(n-1);
Vnom    =   exp(0*1i)*mpc.bus(sb,8);
Vpv     =   mpc.bus(pv,8);
Vmax    =   mpc.bus(:,12);
Vmin    =   mpc.bus(:,13);
bi      =   mpc.branch(:,1);
bj      =   mpc.branch(:,2);



Qlims   =   zeros(length(pv),2);
for it=1:length(pv)
    Qlims(it,1) =   mpc.gen(mpc.gen(:,1)==pv(it),5)-mpc.bus(pv(it),4);
    Qlims(it,2) =   mpc.gen(mpc.gen(:,1)==pv(it),4)-mpc.bus(pv(it),4);
end
Qlims       =   Qlims/mpc.baseMVA;

Qlimsb      =   zeros(2,1);
Qlimsb(1)   =   mpc.gen(mpc.gen(:,1)==sb,5);
Qlimsb(2)   =   mpc.gen(mpc.gen(:,1)==sb,4);

Qlimsb      =   Qlimsb/mpc.baseMVA;
Pmax        =   (mpc.gen(mpc.gen(:,1)==sb,9)-mpc.bus(sb,3))/mpc.baseMVA;


Qss         =   GetQuads(@funpf,nn);
[Qse,Qsi]   =   MakeQs();
Qss         =   ConvertQuads(Qss);
Qse         =   ConvertQuads(Qse);
Qsi         =   ConvertQuads(Qsi);
    


     %   Power Flow Operator %
    function f=funpf(V)
        [Vx,Vy] =   formV(V);
        Vc      =   Vx+1i*Vy;
        S       =   Vc.*conj(Y*Vc);
        f       =   [real(S(nsb))-plims(nsb,1);plims(nsb,2)-real(S(nsb));...
                        imag(S(pq))-qlims(pq,1);qlims(pq,2)-imag(S(pq))];
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
        
            %   Form Voltages and Power Flows   %
        [Vx,Vy] =   formV(V);
        Vc      =   Vx+1i*Vy;
        S       =   Vc.*conj(Y*Vc);  
        
        
            %   PV-bus voltage equality constraints +%
        f       =   abs(Vc(pv)).^2-Vpv.^2;
        
            %   Reactive power limits                %
        f       =   [f;imag(S(pv))-Qlims(:,1);Qlims(:,2)-imag(S(pv))];
        f       =   [f;imag(S(sb))-Qlimsb(1);Qlimsb(2)-imag(S(sb))];
        
            %   Real power limits                %
        f       =   [f;Pmax-real(S(sb))];       
        
            %   Voltage magnitude bounds at PQ buses   %
      %  f       =   [f;abs(Vc(pq)).^2-Vmin(pq).^2];
      %  f       =   [f;Vmax(pq).^2-abs(Vc(pq)).^2];
        
            %   Line flow constraints   %
        f       =   [f;real(Vc(bi).*conj(Vc(bj)))-gam*abs(Vc(bi)).^2;...
                     real(Vc(bi).*conj(Vc(bj)))-gam*abs(Vc(bj)).^2];
       
    end

    %   Constraints on Voltages in Quadratic Form   %
    function [Qseq,Qsiq]=MakeQs
        Qs      =   GetQuads(@func,nn);
        Qseq    =   Qs(:,:,1:length(pv));
        Qsiq    =   Qs(:,:,(length(pv)+1):end);
    end

    %   Function to test correctness of Quadratic Form Constraints   %
    function chk=TestQs(V)
        Ft      =   func(zeros(nn,1));
        Qss     =   cat(3,Qse,Qsi);
        for kk=1:length(Ft)
            Ft(kk)  =   [1;V]'*Qss(:,:,kk)*[1;V];
        end
        chk     =   norm(Ft-func(V));
    end
ff=@fres;
TQs=@TestQs;
end