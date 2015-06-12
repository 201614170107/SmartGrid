function [map,kf,cou]=AddIndices(map,is,js,n,nc)
m       =   length(is);
vc      =   [n^3,n^2,n,1];
cou     =   nc;
for it=1:m
    for jt=it:m
        vecn        =   vecNum([is([it;jt]);js([it;jt])]);
        if(map(vecn)==0)
            cou              =   cou+1;
            map(vecn)        =   cou;
        end
    end
end
function key=vecNum(ks)
ks(ks==0)   =   [];
ks      =   sort(ks(:));
key     =   vc(end-length(ks)+1:end)*ks+1;
end


kf      =   @(inds) map(vecNum(inds));

end
