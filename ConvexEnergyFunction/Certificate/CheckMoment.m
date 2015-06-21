function [opts,res]=CheckMoment(mpc,Scal,gam)
[Qse,Qsi,Qss,Qssp,n]=   MakeMats(mpc,gam,Scal,false);
clqs                =   ComputeGraph([],[Qse;Qsi;Qss]);
[~,kf,ny]           =   FindMap(clqs,n);
QR                  =   QuadRelax(kf);

Mse                 =   QR.transformQuads(Qse,clqs);
Mss                 =   QR.transformQuads(Qss,clqs);
Msi                 =   QR.transformQuads(Qsi,clqs);

opts                =   ones(length(Qsi),1);
[Yo,YI,Yclqs]       =   QR.GetYs(clqs,n);
res                 =   [];

for iter=0:length(Qsi)
    cvx_begin
    variable y(ny,1);
    minimize(0)
    subject to
    if(iter>0)
        QuadPoly.EvalY(Msi{iter},y)==0;
    end
    y(Yo)    ==  1;
    y(YI)    ==  semidefinite(size(YI,1));
    if(Scal.norm==1)
        for it=1:length(Qssp)
            trace(y(YI)'*Qssp{it}.Q)>=0;
        end
    end
    for it=1:length(Qse)
        QuadPoly.EvalY(Mse{it},y)==0;
    end
    for it=1:length(Qsi)
        Yc  =   QuadPoly.EvalY(Msi{it},y);
        Yc  ==  semidefinite(size(Yc,1));
    end
    for it=1:length(Qss)
        Yc  =   QuadPoly.EvalY(Mss{it},y);
        if(~isempty(Qssp))
            Yc  ==  0;
        else
            Yc  ==  semidefinite(size(Yc,1));
        end
    end
    for it=1:length(clqs)
        y(Yclqs{it})    ==  semidefinite(size(Yclqs{it},1));
    end
    cvx_end
    if(iter==0)
        if((strcmp(cvx_status,'Infeasible')||isinf(cvx_optval)))
            opts(:) =   0;
            Qsi     =   [];
            break;
        end
    else
        if(~(strcmp(cvx_status,'Infeasible')||isinf(cvx_optval)))
            opts(iter:end)      =   0;
            Yres                =   y(YI);
            vsol                =   Yres(2:end,1);
            Qe                  =   zeros(length(Qse),1);
            for it=1:length(Qse)
                Qe(it)  =   [1;vsol]'*(Qse{it}.Q)*[1;vsol];
            end
            Qs                  =   zeros(length(Qss),1);
            for it=1:length(Qss)
                Qs(it)  =   [1;vsol]'*(Qss{it}.Q)*[1;vsol];
            end
            Qi      =   zeros(length(Qsi),1);
            for it=1:length(Qsi)
                Qi(it)  =   [1;vsol]'*Qsi{it}.Q*[1;vsol];
            end
            res =   struct('vsol',vsol,'Qs',Qs,'Qe',Qe,'Qi',Qi);
            break;
        end
    end
end
end

