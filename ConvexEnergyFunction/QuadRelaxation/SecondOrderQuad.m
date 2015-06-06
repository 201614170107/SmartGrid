function Y=SecondOrderQuad(Q,clqs,kf,n)
    Ql  =   ConvertStr(Q);
    clq =   FindClq(clqs,Ql,n);
    Y   =   zeros(length(clq)+1)*kf([]);
    for it=1:size(Q,1)
        Y   =   Y+FormMat([0;clq],Q{it,1},Q{it,2},Q{it,3},kf);
    end
end
function Yc=FormMat(clq,qa,qb,qc,kf)
    n   =   length(clq);
    Yc  =   zeros(n)*kf([]);
    for it=1:n
        for jt=it:n
            Yc(it,jt)   =   kf([{num2str(clq(it))},{num2str(clq(jt))},qa,qb]);
            Yc(jt,it)   =   Yc(it,jt);
        end
    end
    Yc  =   Yc*qc;
end
function Ql=ConvertStr(Q)
    Q       =   Q(:,1:2);
    Q       =   Q(:);
    Ql      =   zeros(length(Q),1);
    for it=1:length(Ql)
        if(iscell(Q{it}))
            if(ischar(Q{it}{1}))
                Ql(it)  =   str2double(Q{it}{1});
            end
        end
    end
    Ql  =   unique(Ql);
end
function clo=FindClq(clqs,ii,n)
    ls  =   inf;
    clo =   [];
    ii  =   ii(ii>0);
    for it=1:length(clqs)
        clq =   clqs{it};
        v   =   sparse(clq,1,1,n,1);
        if(all(v(ii)>0))
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