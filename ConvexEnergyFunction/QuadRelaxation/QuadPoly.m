classdef QuadPoly
    properties
        is;
        js;
        vs;
        n;
        Q;
    end
    methods
        function obj = QuadPoly(Q,vname,nn)
            Q(abs(Q)<1e-5*max(abs(Q(:))))   =   0;
            obj.Q      =   Q;
            [i,j,v]    =   find(Q);
            ii         =   i>j;
            ie         =   i==j;
            obj.is     =   vname([i(ii);i(ie)]);
            obj.js     =   vname([j(ii);j(ie)]);
            obj.vs     =   [2*v(ii);v(ie)];
            obj.n      =   nn;
        end
        function Q  =   GetQ(obj)
            Q   =   obj.Q;
        end
        function Ms =   FormMat(obj,clqs,GetY)
            if(iscell(clqs))
                inds   =   obj.FindClq(clqs);
                inds   =   [0;inds];
            else
                inds   =    clqs;
            end
            nz     =   length(inds);
            Ms     =   cell(length(obj.is),1);
            for i=1:length(obj.is)
                Mcv =   zeros(nz);
                Mcc =   zeros(nz);                    
                for it=1:nz
                    for jt=it:nz
                        [Mcv(it,jt),Mcc(it,jt)]  =    GetY([{{inds(it)}},...
                                                             {{inds(jt)}},...
                                                             obj.is(i),...
                                                             obj.js(i)]);
                        Mcv(jt,it)  =   Mcv(it,jt);
                        Mcc(jt,it)  =   Mcc(it,jt);
                    end
                end
                Mcv     =   Mcv*obj.vs(i);
                Ms{i}   =   {Mcv,Mcc};
            end
        end
        
        
        function clo   =   FindClq(obj,clqs)
            ind     =   obj.GetInds();
            ls      =   inf;
            clo     =   [];
            for it=1:length(clqs)
                clq =   clqs{it};
                cl  =   sparse(clq,1,1,obj.n,1)>0;
                if(all(cl(ind)))
                    if(length(clq)<ls)
                        ls  =   length(clq);
                        clo =   clq;
                    end
                end
            end
            if(isempty(clo))
                disp('trouble');
            end
        end
        
        function ind=GetInds(obj)
            i   =   arrayfun(@(anyx)toind(anyx{1}),obj.is);
            j   =   arrayfun(@(anyx)toind(anyx{1}),obj.js);                  
            ind     =   unique([i;j]);
            ind(ind==0) =   [];
        end
        
        function res=Eval(obj,y)
            res=[1;y]'*(obj.Q)*[1;y];
        end
    end
    methods(Static)
        function Ym=EvalY(Ms,y)
            Ym  =   Ms{1}{1}.*y(Ms{1}{2});
            for it=2:length(Ms)
                Ym  =   Ym+Ms{it}{1}.*y(Ms{it}{2});
            end
        end
    end
end

function res=toind(anyx)
res=0;
if(iscell(anyx))
    res=anyx{1};
end
end