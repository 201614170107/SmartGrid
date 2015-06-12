classdef QuadMatPoly
    properties
        M;
        zi;
        Qp;
    end
    methods
         function obj = QuadMatPoly(Ms,zz,Q)
            obj.Qp  =   Q;
            obj.M   =   Ms;
            obj.zi  =   zz;
         end
         function res = FormPol(obj,GetY)
             is         =   obj.zi;
             res        =   obj.Qp.FormMat(is,GetY);
             res        =   {is,res,obj.M};
         end
    end
    methods (Static)
        function pol=EvalPoly(res,y,n)
            pol =   zeros(n)*y(1);
            for it=1:length(res)
                is  =   res{it}{1};
                Mc  =   res{it}{3};                                
                resc=   res{it}{2};
                for jt=1:length(resc)
                    Msv =   resc{jt}{1};
                    Msi =   resc{jt}{2};
                    pol(is,is)  =   pol(is,is)+Mc*(Msv.*y(Msi));
                end
            end
        end
    end
end