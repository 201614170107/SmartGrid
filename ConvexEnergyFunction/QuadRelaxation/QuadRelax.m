classdef QuadRelax
    properties
        kfun;
    end
    methods
        function obj=QuadRelax(kf)
            obj.kfun  =   kf;
        end
        function [v,i]=GetYkf(obj,inds)
            kf  =   obj.kfun;
            if(isempty(inds))
                v   =   1;
                i   =   1;
            else
                m=1;
                indv=zeros(length(inds),1);
                for it=1:length(indv)
                    if(iscell(inds{it}))
                        indv(it) =   inds{it}{1};
                    else
                        m       =   m*inds{it};
                        indv(it) =   0;
                    end
                end
                v   =   m;
                i   =   kf(indv);
                if(i==0)
                    disp('trouble');
                end
            end
        end
        
        function Ms=transformQuads(obj,Qs,clqs)
            GY  =   @(anyinds) obj.GetYkf(anyinds);
            Ms  =   cell(length(Qs),1);
            for it=1:length(Qs)
                Ms{it}  =   Qs{it}.FormMat(clqs,GY);
            end
        end
        
        function Ms=transformQuadMats(obj,Mats)
            GY  =   @(anyinds) obj.GetYkf(anyinds);
            Ms  =   cell(length(Mats),1);
            for it=1:length(Ms)
                Ms{it}  =   Mats{it}.FormPol(GY);
            end
        end
        
        function [is,js]=FormPairs(~,ii)
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
        
        function [Yo,YI,Yclqs]=GetYs(obj,clqs,n)
            Yo                  =   obj.kfun([0,0,0,0]);
            YI                  =   obj.FormMkf(n);
            Yclqs               =   cell(length(clqs),1);
            for it=1:length(clqs)
            [is,js]         =   obj.FormPairs(clqs{it});
            Yclqs{it}       =   obj.FormMatkf(is,js);
            end
        end
        
        function Yc=FormMatkf(obj,is,js)
            kf  =   obj.kfun;
            Yc  =   zeros(length(is));
            for it=1:length(is)
                for jt=it:length(js)
                    Yc(it,jt)   =   kf([is(it),js(it),is(jt),js(jt)]);
                    Yc(jt,it)   =   Yc(it,jt);
                end
            end
        end
        
        function Yc=FormMkf(obj,n)
            kf  =   obj.kfun;
            Yc  =   zeros(n+1);
            for it=1:(n+1)
                for jt=it:(n+1)
                    Yc(it,jt)   =   kf([it-1,jt-1,0,0]);
                    Yc(jt,it)   =   Yc(it,jt);
                end
            end
        end
        
    end
end