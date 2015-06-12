function [Qse,Qbr]=FormQuads(mpc,gam)
pv  =   find(mpc.bus(:,2)==2);
pq  =   find(mpc.bus(:,2)==1);
sb  =   find(mpc.bus(:,2)==3);
nsb =   find(~(mpc.bus(:,2)==3));
nbus=   size(mpc.bus,1);
bi  =   mpc.branch(:,1);
bj  =   mpc.branch(:,2);
Vpv =   mpc.bus(pv,8);
npv =   length(pv);
Qs  =   GetQuads(@QuadFuns,2*(nbus-1));
m   =   size(Qs,3);
Qse =   cell(npv,1);
Qbr =   cell(m-npv,1);

v   =   (1:2*(nbus-1))'+(nbus-1+length(pq));
v   =   arrayfun(@(anyx) {num2str(anyx)},v);
v   =   [1;v];

n   =   3*(nbus-1)+length(pq);
for it=1:(m-npv)
    Qbr{it}         =   QuadPoly(Qs(:,:,it),v,n);
end
for it=(m-npv+1):m
    Qse{it-(m-npv)} =   QuadPoly(Qs(:,:,it),v,n);
end

    function V=formV(v)
        V           =   zeros(nbus,1);
        V(nsb)      =   v(1:(nbus-1))+1i*v(nbus:end);
        V(sb)       =   mpc.bus(sb,8);       
    end
    function f=QuadFuns(v)
        Vc   =  formV(v);
        f    =  [real(Vc(bi).*conj(Vc(bj)))-gam*abs(Vc(bi)).^2;...
                 real(Vc(bi).*conj(Vc(bj)))-gam*abs(Vc(bj)).^2;...
                 abs(Vc(pq)).^2-.6.^2;...
                 abs(Vc(pv)).^2-Vpv.^2];
    end
end

