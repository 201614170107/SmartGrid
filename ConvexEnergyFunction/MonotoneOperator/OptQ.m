function Q=OptQ(anyM,P)

Sym =   @(anyx) anyx+anyx';
n=size(anyM,1);
        cvx_begin
            variable Q(n,n);
            maximize(lambda_min(Sym(P*anyM*Q)))
            subject to
            norm(Q)<=1
        cvx_end
      %  Q=diag(Q);
    end
