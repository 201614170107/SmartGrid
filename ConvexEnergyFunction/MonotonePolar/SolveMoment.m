function opt=SolveMoment(Mats,mpc,gam,opt)
nbus            =   size(mpc.bus,1);
npq             =   sum(mpc.bus(:,2)==1);
n               =   3*(nbus-1)+npq;
nz              =   (nbus-1)+npq;
[Qse,Qsi]       =   FormQuads(mpc,gam,opt);
clqs            =   ComputeGraph(Mats,[Qse;Qsi]);
[~,kf,ny]       =   FindMap(clqs,n);
QR              =   QuadRelax(kf);
Mse             =   QR.transformQuads(Qse,clqs);
Msi             =   QR.transformQuads(Qsi,clqs);
[Yo,YI,Yclqs]   =   QR.GetYs(clqs,n);
Mats            =   QR.transformQuadMats(Mats);

z               =   (1:nz)';
cvx_begin
variable y(ny,1);
pol         =   QuadMatPoly.EvalPoly(Mats,y,nz);
minimize(0)
subject to
-(pol+pol')         ==  semidefinite(nz);
y(Yo)               ==  1;
y(YI)               ==  semidefinite(size(YI,1));
normz               =   0;
for it=1:nz
    ii      =   z(it);
    normz   =   normz+y(kf([ii,ii,0,0]));
end
normz   ==  1;
for it=1:length(Qse)
    QuadPoly.EvalY(Mse{it},y)==0;
end
for it=1:length(Qsi)
    Yc  =   QuadPoly.EvalY(Msi{it},y);
    Yc  ==  semidefinite(size(Yc,1));
end
for it=1:length(clqs)
    y(Yclqs{it})  ==  semidefinite(size(Yclqs{it},1));
end
cvx_end
opt=cvx_optval;

end
