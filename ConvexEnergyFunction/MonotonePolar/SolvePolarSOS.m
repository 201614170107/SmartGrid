function sol=SolvePolarSOS(Mv,Mbr,mpc,del,vscl,vr)
define_constants;
pq      =   mpc.bus(:,BUS_TYPE)==1;
pv      =   mpc.bus(:,BUS_TYPE)==2;
sb      =   mpc.bus(:,BUS_TYPE)==3;
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
v       =   msspoly('v',length(Mv));
s       =   msspoly('s',1);
c       =   msspoly('c',1);



prog    =   prog.withIndeterminate([z;v;s;c]);
V       =   zeros(nbus,1)*v(1);
V(sb)   =   mpc.bus(sb,VM);
V(pv)   =   mpc.bus(pv,VM);
V(pq)   =   v;
[prog,W]=   prog.newFree(n);
pol     =   [];
Sym     =   @(anyx,indx) (diag(W(indx))*anyx)+...
    (diag(W(indx))*anyx)';
Xs      =   cell(nbus,1);
for it=1:nbus
    Xs{it}  =   zeros(n)*W(1);
end

for it=1:length(Mbr)
    Ms          =   Mbr{it}{1};
    Mc          =   Mbr{it}{2};
    ind         =   Mbr{it}{3};
    [prog,X]    =   ApplyRelax(prog,Sym(Ms,ind),Sym(Mc,ind),[z(ind);s;c],del);
    if(isempty(pol))
        pol     =   z(ind)'*(X*V(bi(it))*V(bj(it)))*z(ind);
    else
        pol     =   pol+z(ind)'*(X*V(bi(it))*V(bj(it)))*z(ind);
    end
end

for it=1:length(Mv)
    Mx  =   Mv{it}{1};
    ind =   Mv{it}{2};
    ii  =   Mv{it}{3};
    pol =   pol+V(ii)^2*z(ind)'*Sym(Mx,ind)*z(ind);
end

mbr                 =   pq(bi)|pq(bj);
Qsi                 =   MakeQuad([v-Vmin(pq);Vmax(pq)-v]);
%Qsi                 =   [Qsi;(V(bi(mbr))-vr*V(bj(mbr)));...
 %                           (V(bj(mbr))-vr*V(bi(mbr)));...
  %                          (V(bi(mbr))-vr*V(bj(mbr))).*(V(bj(mbr))-vr*V(bi(mbr)));...
  %                          (V(bi(mbr))-vr*V(bj(mbr))).^2;...
   %                         (V(bj(mbr))-vr*V(bi(mbr))).^2];
   %...
                   %         V(bi(mbr))-vr*V(bj(mbr));...
                    %        V(bj(mbr))-vr*V(bi(mbr));
Msi                 =   Qsi;
for it=1:length(Qsi)
[prog,tmp]      =   prog.newPSD(1+length([z]));
Msi(it)         =   [1;z]'*tmp*[1;z];      
end
[prog,Mse]          =   prog.newFreePoly(monomials([z;v],0:2),1);
prog                =   prog.withSOS(pol-Msi'*Qsi-Mse*(z'*z-1)-1);

options             =   spot_sdp_default_options();
options.verbose     =   1;
sol                 =   prog.minimize(0,@spot_frlib,options);
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
