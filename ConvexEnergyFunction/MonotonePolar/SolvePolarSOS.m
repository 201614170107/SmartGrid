function sol=SolvePolarSOS(Mv,Mbr,mpc,vscl,gam)
define_constants;
pq      =   mpc.bus(:,BUS_TYPE)==1;
pv      =   mpc.bus(:,BUS_TYPE)==2;
sb      =   mpc.bus(:,BUS_TYPE)==3;
nsb     =   ~sb;
Vmin    =   mpc.bus(:,13);
Vmax    =   mpc.bus(:,12);
Vmean   =   (Vmin+Vmax)/2;
Vmin    =   Vmean+(Vmin-Vmean)*vscl;
Vmax    =   Vmean+(Vmax-Vmean)*vscl;

npq     =   sum(double(pq));
bi      =   mpc.branch(:,1);
bj      =   mpc.branch(:,2);
prog    =   spotsosprog;
nbus    =   size(mpc.bus,1);
n       =   nbus+npq-1;
z       =   msspoly('z',n);
vxv     =   msspoly('vx',(nbus-1));
vyv     =   msspoly('vy',(nbus-1));

prog    =   prog.withIndeterminate([z;vxv;vyv]);
V       =   zeros(nbus,2)*vxv(1);
V(find(sb),1) =   mpc.bus(sb,VM);
V(find(sb),2) =   0;
V(find(nsb),1)=   vxv;
V(find(nsb),2)=   vyv;
%[prog,W]=   prog.newFree(1);
pol     =   [];
Sym     =   @(anyx,indx) (anyx+anyx');

for it=1:length(Mbr)
    Ms          =   Mbr{it}{1};
    Mc          =   Mbr{it}{2};
    ind         =   Mbr{it}{3};
    ic          =   bi(it);
    jc          =   bj(it);
    pols        =   z(ind)'*(Sym(Ms,ind)*(V(jc,1)*V(ic,2)-V(ic,1)*V(jc,2)))*z(ind);
    polc        =   z(ind)'*(Sym(Mc,ind)*(V(ic,1)*V(jc,1)+V(ic,2)*V(jc,2)))*z(ind);
    if(isempty(pol))
        pol     =   pols+polc;
    else
        pol     =   pol+pols+polc;
    end
end

for it=1:length(Mv)
    Mx  =   Mv{it}{1};
    ind =   Mv{it}{2};
    ii  =   Mv{it}{3};
    pol =   pol+(V(ii,1).^2+V(ii,2).^2)*z(ind)'*Sym(Mx,ind)*z(ind);
end

Qse     =   [V(pv,1).^2+V(pv,2).^2-mpc.bus(pv,VM).^2;z'*z-1];

Qsi     =   V(bi,1).*V(bj,1)+V(bi,2).*V(bj,2)-gam*(V(bi,1).^2+V(bi,2).^2);
Qsi     =   [Qsi;V(bi,1).*V(bj,1)+V(bi,2).*V(bj,2)-gam*(V(bj,1).^2+V(bj,2).^2)];

[prog,Msi]     =   MakeMult(prog,length(Qsi),[z;vxv;vyv],1);
[prog,Mse]     =   MakeMult(prog,length(Qse),[z;vxv;vyv],0);
prog    =   prog.withSOS(pol-Msi'*Qsi-Mse'*Qse);

options             =   spot_sdp_default_options();
options.verbose     =   1;
sol                 =   prog.minimize(0,@spot_frlib,options);
end

function [prog,Mult]=MakeMult(prog,li,var,pos)
    Mult    =   [];
    for i=1:li
        if(pos)
            [prog,tmp]  =   prog.newPSD(length(var)+1);
            Mult        =   [Mult;[1;var]'*tmp*[1;var]];
        else
           [prog,tmp]   =   prog.newSym(length(var)+1);
            Mult        =   [Mult;[1;var]'*tmp*[1;var]];
        end        
    end       
end

function [prog,X]=ApplyRelax(prog,Ms,Mc,vars,del)
n           =   size(Ms,1);
z           =   vars(1:n);
s           =   vars(n+1);
c           =   vars(n+2);
[prog,X]    =   prog.newSym(n);
pol         =   z'*(Ms*s+Mc*c-X)*z;


zz          =   vars;
Qsi         =   MakeQuad([c-cos(del);1-c;sin(del)-s;s+sin(del)]);
Msi         =   Qsi;
for it=1:length(Qsi)
[prog,tmp]  =   prog.newPSD(1+length(z));
Msi(it)     =   [1;z]'*tmp*[1;z];
end
Qse         =   [s^2+c^2-1;z'*z-1];
[prog,Mse]  =   prog.newFreePoly(monomials(zz,0:2),length(Qse));

prog        =   prog.withSOS(pol-Qsi'*Msi-Qse'*Mse);

end

function Qs=MakeQuad(Qs)
    Qs  =   [1;Qs];
    Qs  =   Qs*Qs';
    Qs  =   Qs(triu(true(size(Qs,1))));
    Qs  =   Qs(2:end);
end

function xys=formxys(del,levs)
dels    =   linspace(-del,del,2*levs+1)';
M       =    [0,1,cos(del);...
    sin(dels),cos(dels),ones(size(dels));...
    0,1,cos(del)];
xys =   zeros(size(M,1)-1,2);
for it=1:size(M,1)-1
    xys(it,:)   =   (M(it:(it+1),1:2)\M(it:(it+1),3))';
end
end

function [prog,M]=newASym(prog,n)
if(n>1)
    [prog,M]    =   prog.newFree(n*n);
    M           =   reshape(M,[n,n]);
    prog        =   prog.withEqs(M+M');
else
    M=0;
end
end
