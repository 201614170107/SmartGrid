function [Qse,Qsi,Qss,Qssp,nn,ff]=MakeMats(mpc,gam,Scal,lflim)
Y       =   makeYbus(mpc);
sb      =   find(mpc.bus(:,2)==3);
n       =   size(Y,1);
nsb     =   setdiff(1:n,sb);
pv      =   find(mpc.bus(:,2)==2);
pq      =   setdiff(nsb,pv);
nn      =   2*(n-1);
if(Scal.norm==1)
    nn  =   nn+length(nsb)+length(pq);
end
mpcres  =   runpf(mpc);
Vnom    =   exp(0*1i)*mpcres.bus(sb,8);
Vpv     =   mpcres.bus(pv,8);
Vmax    =   mpc.bus(:,12);
Vmin    =   mpc.bus(:,13);
Y       =   makeYbus(mpc);
bi      =   mpc.branch(:,1);
bj      =   mpc.branch(:,2);
llim    =   [];
if(lflim)
    y       =   Y(sub2ind(size(Y),bi,bj));
    llim    =   [bi,bj,(mpc.branch(:,6)./abs(y))./mpc.baseMVA];
end



Qlims   =   zeros(length(pv),2);
for it=1:length(pv)
    Qlims(it,1) =   mpc.gen(mpc.gen(:,1)==pv(it),5)-mpc.bus(pv(it),4);
    Qlims(it,2) =   mpc.gen(mpc.gen(:,1)==pv(it),4)-mpc.bus(pv(it),4);
end
Qlims       =   Qlims/mpc.baseMVA;

Qlimsb      =   [mpc.gen(mpc.gen(:,1)==sb,5)-mpc.bus(sb,4);...
                 mpc.gen(mpc.gen(:,1)==sb,4)-mpc.bus(sb,4)]/mpc.baseMVA;
Plimsb      =   [mpc.gen(mpc.gen(:,1)==sb,10)-mpc.bus(sb,4);...
                 mpc.gen(mpc.gen(:,1)==sb,9)-mpc.bus(sb,4)]/mpc.baseMVA;


Qss         =   GetQuads(@funpf,nn);
[Qse,Qsi]   =   MakeQs();
Qse         =   ConvertQuads(Qse);
Qss         =   ConvertQuads(Qss);
Qsi         =   ConvertQuads(Qsi);
if(Scal.norm==1)
    Qssp     =   GetQuads(@sConst,nn);
    Qssp     =   ConvertQuads(Qssp);
else
    Qssp     =  [];
end


     %   Power Flow Operator %
    function f=funpf(V)
        [Vx,Vy,s] =   formV(V);
        Vc      =   Vx+1i*Vy;
        S       =   Vc.*conj(Y*Vc);
        if(Scal.norm==0)
            plims   =   Scal.Getplims;
            qlims   =   Scal.Getqlims;
            f       =   [real(S(nsb))-plims(:,1);plims(:,2)-real(S(nsb));...
                        imag(S(pq))-qlims(:,1);qlims(:,2)-imag(S(pq))];
        else
           f       =   [real(S(nsb))-s(1:(n-1));imag(S(pq))-s(n:end)];
        end
    end

    %   Helper Function - Add Slack Bus Voltage Reference %
    function [Vx,Vy,s]=formV(V)
        Vx      =   ones(n,1)*real(Vnom);
        Vy      =   ones(n,1)*imag(Vnom);
        Vx(nsb) =   V(1:(n-1))+real(Vnom);
        Vy(nsb) =   V(n:(2*(n-1)))+imag(Vnom);
        if(Scal.norm==1)
            s       =   V(2*n-1:end);
        else
            s       =   [];
        end
    end

    %   Helper Function - Add Slack Bus Voltage Reference %
    function res=sConst(V)
        [~,~,s] =   formV(V);
        res     =   Scal.MakeQuad(s);
    end

    %   Constraints on Voltages %
    function f=func(V)
        
            %   Form Voltages and Power Flows   %
        [Vx,Vy] =   formV(V);
        Vc      =   Vx+1i*Vy;
        S       =   Vc.*conj(Y*Vc);  
        
        
            %   PV-bus voltage equality constraints +%
        f       =   abs(Vc(pv)).^2-Vpv.^2;
        
        %{
            %   Reactive power limits                %
        f       =   [f;Qlims(:,2)-imag(S(pv));imag(S(pv))-Qlims(:,1)];
        f       =   [f;Qlimsb(2)-imag(S(sb));imag(S(sb))-Qlimsb(1)];
        
            %   Real power limits                %
        f       =   [f;Plimsb(2)-real(S(sb));real(S(sb))-Plimsb(1)];   
        %}
        %   Voltage magnitude bounds at PQ buses   %
        f       =   [f;abs(Vc(pq)).^2-Vmin(pq).^2];
        f       =   [f;Vmax(pq).^2-abs(Vc(pq)).^2];
        
            %   Line flow constraints   %
        f       =   [f;gam.^2-abs(Vc(bi)-Vc(bj)).^2];
        if(~isempty(llim))
            f   =   [f;-abs(V(llim(:,1))-V(llim(:,2))).^2+llim(:,3).^2];
        end
       %}
    end

    %   Constraints on Voltages in Quadratic Form   %
    function [Qseq,Qsiq]=MakeQs
        Qs      =   GetQuads(@func,nn);
        Qseq    =   Qs(:,:,1:length(pv));
        Qsiq    =   Qs(:,:,(length(pv)+1):end);
    end
%{
    %   Function to test correctness of Quadratic Form Constraints   %
    function chk=TestQs(V)
        Ft      =   func(zeros(nn,1));
        Qss     =   cat(3,Qse,Qsi);
        for kk=1:length(Ft)
            Ft(kk)  =   [1;V]'*Qss(:,:,kk)*[1;V];
        end
        chk     =   norm(Ft-func(V));
    end
%}
ff=@fres;
%TQs=@TestQs;
end