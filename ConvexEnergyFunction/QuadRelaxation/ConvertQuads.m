function Qsn=ConvertQuads(Qs)
    n   =   size(Qs,1);
    m   =   size(Qs,3);
    x   =   arrayfun(@(anyx) {{anyx}},(1:(n-1))');
    x   =   [1;x];
    Qsn =   cell(m,1);
    for it=1:m
        Q           =   Sparsify(Qs(:,:,it));
        Qsn{it}     =   QuadPoly(Q,x,n-1);
    end
end
function Q=Sparsify(Q)
    Q(abs(Q)<1e-5*max(max(abs(Q))))  =   0;
    Q                           =   sparse(Q);
end