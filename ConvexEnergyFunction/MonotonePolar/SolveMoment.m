function opt=SolveMoment(Mats,mpc,gam)
nbus            =   size(mpc.bus,1);
npq             =   sum(mpc.bus(:,2)==1);
n               =   3*(nbus-1)+npq;
nz              =   (nbus-1)+npq;
[Qse,Qsi]       =   FormQuads(mpc,gam);
clqs            =   ComputeGraph(mpc,Mats,Qse,Qsi);
[~,kf,ny]       =   FindMap(clqs,n);
QR              =   QuadRelax(kf);
Mse             =   QR.transformQuads(Qse);
Msi             =   QR.transformQuads(Qsi);
z               =   (1:nz)';
cvx_begin
variable y(ny,1);
GetY        =   @(anyind) GetYkf(anyind,y,kf);
pol         =   zeros(nz)*y(1);
for it=1:length(Mats)
    pol         =   Mats{it}.UpdatePol(pol,GetY);
end
minimize(0)
subject to
-(pol+pol')         ==  semidefinite(size(pol,1));
y(kf([0,0,0,0]))    ==  1;
normz=0;
for it=1:nz
    ii      =   z(it);
    normz   =   normz+y(kf([ii,ii,0,0]));
end
normz   ==  1;
y(FormMkf(n,kf))==semidefinite(n+1);
for it=1:length(Qse)
    Qse{it}.FormMat(clqs,GetY)==0;
end
for it=1:length(Qsi)
    Yc  =   Qsi{it}.FormMat(clqs,GetY);
    Yc  ==  semidefinite(size(Yc,1));
end
for it=1:length(clqs)
    [is,js] =   FormPairs(clqs{it});
    Yc  =   y(FormMatkf(is,js,kf));
    Yc  ==  semidefinite(size(Yc,1));
end
cvx_end
opt=cvx_optval;

end

function res=GetYkf(inds,y,kf)
if(isempty(inds))
    res=y(1);
else
    m=1;
    indv=zeros(length(inds),1);
    for it=1:length(indv)
        if(ischar(inds{it}))
            indv(it) =   str2double(inds{it});
        else
            m       =   m*inds{it};
            indv(it) =   0;
        end
    end
    res=m*y(kf(indv));
    if(isempty(res))
        disp('trouble');
    end
end
end

function [is,js]=FormPairs(ii)
    ii  =   [0;ii];
    n   =   length(ii);
    is  =   zeros(n*(n+1)/2,1);
    js  =   is;
    cc  =   0;
    for it=1:n
        for jt=it:n
            cc      =   cc+1;
            is(cc)  =   ii(it);
            js(cc)  =   ii(jt);
        end
    end
end

function Yc=FormMatkf(is,js,kf)
Yc  =   zeros(length(is));
for it=1:length(is)
    for jt=it:length(js)
        Yc(it,jt)   =   kf([is(it),js(it),is(jt),js(jt)]);
        Yc(jt,it)   =   Yc(it,jt);
    end
end
end

function Yc=FormMkf(n,kf)
Yc  =   zeros(n+1);
for it=1:(n+1)
    for jt=it:(n+1)
        Yc(it,jt)   =   kf([it-1,jt-1,0,0]);
        Yc(jt,it)   =   Yc(it,jt);
    end
end
end
