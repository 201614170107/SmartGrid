classdef Sset
    properties
        Sigma;
        snom;
        rad;
        plims;        
        qlims;
        norm;
    end
    methods
        function obj=Sset(norm,rad,va,vb)
            if(norm==0)
                obj.norm    =   norm;                
                obj.Sigma   =   [];
                obj.snom    =   [];
                obj.rad     =   rad;
                obj.plims   =   va;
                obj.qlims   =   vb;
            end
            if(norm==1)
                obj.plims   =   [];
                obj.qlims   =   [];
                obj.norm    =   norm;
                obj.Sigma   =   va;
                obj.rad     =   rad;
                obj.snom    =   vb;
            end
        end
    function Q=MakeQuad(obj,V)
        bd  =   [obj.plims;obj.qlims];
        if(obj.norm==0)
            Q   =   Quadbd(V,bd);
        else
            Q   =   obj.rad.^2-sum(((obj.Sigma)*(V-obj.snom)).^2);
        end
    end
    function bd=Getplims(obj)
            bd  =   obj.plims;
            bd  =    mean(bd,2)*ones(1,2)+obj.rad*[(-mean(bd,2)+bd(:,1)),...
                                                    (-mean(bd,2)+bd(:,2))];
    
    end
    function bd=Getqlims(obj)
            bd  =   obj.qlims;
            bd  =    mean(bd,2)*ones(1,2)+obj.rad*[(-mean(bd,2)+bd(:,1)),...
                                                    (-mean(bd,2)+bd(:,2))];
    
    end
    end
    
end
function f=Quadbd(x,bds)
    f   =   (x-bds(:,1)).*(bds(:,2)-x);    
end