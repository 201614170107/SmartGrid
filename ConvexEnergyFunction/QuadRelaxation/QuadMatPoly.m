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
         function pol = UpdatePol(obj,pol,GetY)
             is         =   obj.zi;
             Ms         =   obj.Qp.FormMat(is,GetY);
             pol(is,is) =   pol(is,is)+(obj.M)*Ms;
         end
    end
end