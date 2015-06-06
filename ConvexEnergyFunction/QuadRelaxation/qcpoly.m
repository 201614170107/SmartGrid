classdef qcpoly
    properties
        n; %Number of indeterminates
        list; %Cell array with list of monomials        
        sz; %Size of qcpoly object
    end
    methods
        function obj=qcpoly(nc)
            if(nargin>0)
                if(isinteger(nc)&&(nc>0))
                    obj.n=nc;
                    list=cell(nc,1);
                    for it=1:nc
                        list{it}={it};
                    end
                else
                    error('Value must be positive integer');
                end
            end
        end
        function obj=Mult(obj1,obj2)
            if(size(obj1,2)==size(obj2,1))
                for it=1:
            end
        end
        function siz=size(obj)
            siz=sz;
        end
        function objn=subsref(A,s)
            objn=qcpoly(n);
            
         end
    end
end