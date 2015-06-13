function opts=CheckMoment(mpc,plims,qlims,gam)
[Qse,Qsi,Qss,n]     =   MakeMats(mpc,gam,plims,qlims);
clqs                =   ComputeGraph([],[Qse;Qsi;Qss]);
[~,kf,ny]           =   FindMap(clqs,n);
QR                  =   QuadRelax(kf);

Mse                 =   QR.transformQuads(Qse,clqs);
Mss                 =   QR.transformQuads(Qss,clqs);
Msi                 =   QR.transformQuads(Qsi,clqs);
opts                =   ones(length(Qsi),1);
[Yo,YI,Yclqs]       =   QR.GetYs(clqs,n);

for iter=1:length(Qsi)
cvx_begin
variable y(ny,1);
minimize(0)
subject to
trace(Qsi{iter}.GetQ()'*y(YI))==0;
y(Yo)    ==  1;
y(YI)    ==  semidefinite(size(YI,1));
for it=1:length(Qse)
    QuadPoly.EvalY(Mse{it},y)==0;
end
for it=1:length(Qsi)
    Yc  =   QuadPoly.EvalY(Msi{it},y);
    Yc  ==  semidefinite(size(Yc,1));
end
for it=1:length(Qss)
    Yc  =   QuadPoly.EvalY(Mss{it},y);    
    Yc  ==  semidefinite(size(Yc,1));
end
for it=1:length(clqs)
    y(Yclqs{it})    ==  semidefinite(size(Yclqs{it},1));
end
cvx_end
if(~strcmp(cvx_status,'Infeasible'))
    opts(:) =   0;
    Qsi=[];
    break;
end
end
%}
end

