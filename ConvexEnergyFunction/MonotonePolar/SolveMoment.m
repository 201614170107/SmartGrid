function opt=SolveMoment(Mv,Mbr,mpc,gam,del)
nbus            =   size(mpc.bus,1);
pq              =   find(mpc.bus(:,2)==1);
pq              =   pq(:);
npq             =   length(pq);
[Qse,Qsi,V]     =   FormQuads(mpc,gam,del);
z               =   arrayfun(@(anyx) {num2str(anyx)},(1:(nbus-1+npq))');
n               =   3*(nbus-1)+npq;
clqs            =   ComputeGraph(mpc,Mv,Mbr,Qse,Qsi);
[~,kf,cou]      =   FindMap(clqs,n);
Mats            =   ConvertMats(Mv,Mbr,V,z);

cvx_begin
variable y(cou,1);
GetY        =   @(anyind) GetYkf(anyind,y,kf);
pol         =   FormPoly(Mats,length(z),GetY);
minimize(lambda_max(pol+pol'))
subject to
y(kf([0,0,0,0]))==1;
normz=0;
for it=1:length(z)
    ii=str2double(z{it});
    normz=normz+y(kf([ii,ii,0,0]));
end
normz   ==  1;
y(FormMkf(n,kf))==semidefinite(n+1);
for it=1:length(Qse)
    SecondOrderQuad(Qse{it},clqs,GetY,n)==0;
end
for it=1:length(Qsi)
    Yc  =   SecondOrderQuad(Qsi{it},clqs,GetY,n);
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

function Mats=ConvertMats(Mv,Mbr,V,z)
    Mats    =   cell(2*length(Mv)+4*length(Mbr),1);
    cou     =   0;
    for it=1:length(Mv)
        M           =   Mv{it}{1};
        ind         =   Mv{it}{2};
        ii          =   Mv{it}{3};
        Mats{cou+1} =   {M,z(ind),V(ii,1),V(ii,1)};
        Mats{cou+2} =   {M,z(ind),V(ii,2),V(ii,2)};
        cou         =   cou+2;
    end
    for it=1:length(Mbr)
        Ms          =   Mbr{it}{1};
        Mc          =   Mbr{it}{2};
        ind         =   Mbr{it}{3};
        ij          =   Mbr{it}{4};
        Mats{cou+1} =   {Mc,z(ind),V(ij(1),1),V(ij(2),1)};
        Mats{cou+2} =   {Mc,z(ind),V(ij(1),2),V(ij(2),2)};
        Mats{cou+3} =   {Ms,z(ind),V(ij(1),2),V(ij(2),1)};
        Mats{cou+4} =   {-Ms,z(ind),V(ij(1),1),V(ij(2),2)};
        cou         =   cou+4;
    end
end

