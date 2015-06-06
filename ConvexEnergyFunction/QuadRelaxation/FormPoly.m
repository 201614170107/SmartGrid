function pol=FormPoly(Mats,n,kf)
    pol =   zeros(n,n)*kf([]);
    for it=1:length(Mats)
        M   =   Mats{it}{1};
        z   =   Mats{it}{2};
        va  =   Mats{it}{3};
        vb  =   Mats{it}{4};
        yc  =   formy(z,va,vb,kf);
        ind =   zeros(length(z),1);
        for jt=1:length(ind)
            ind(jt) =   str2double(z{jt});
        end
        if(length(ind)>1)
            pol(ind,ind) =   pol(ind,ind)+M*yc;        
        else
            pol(ind,ind) =   pol(ind,ind)+M*yc;        
        end
    end
end
function M=formy(z,va,vb,kf)
    n   =   length(z);
    M   =   zeros(n)*kf([]);
    for it=1:n
        for jt=it:n
            M(it,jt)    =   kf([z(it),z(jt),va,vb]);
            M(jt,it)    =   M(it,jt);
        end
    end
end