function [M,cou,kf]=ComputeIndices(is,js,n)
map     =   java.util.HashMap;
m       =   length(is);
cou     =   0;
M       =   zeros(m);
for it=1:m
    for jt=it:m
        vecn        =   vecNum(n,[is([it;jt]);js([it;jt])]);
        if(~map.containsKey(vecn))
            cou              =   cou+1;
            map.put(vecn,cou);
        end
        M(it,jt)    =   map.get(vecn);
        M(jt,it)    =   map.get(vecn);
    end
end


kf      =   @(inds) map.get(vecNum(n,inds));
end

function key=vecNum(n,ks)
vec=sparse(ks+1,1,1,n+1,1);
vec(1)=[];
key=char(vec+32);
end