function [map,kf,cou]=FindMap(clqs,n)
map     =   sparse([],[],[],(n+1)^4,1);
cou     =   0;
for it=1:length(clqs)
    [is,js]         =   FormPairs(clqs{it});
    [map,kf,cou]    =   AddIndices(map,is,js,n,cou);
end
is=(0:n)';
js=zeros(n+1,1);
[map,kf,cou]    =   AddIndices(map,is,js,n,cou);
%M               =   FormTable(kf,n);
end


function [is,js]=FormPairs(ii)
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