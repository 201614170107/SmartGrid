function P=OptP(anyM,Q)

Sym =   @(anyx) anyx+anyx';
n=size(anyM,1);
        cvx_begin
            variable P(n,n);
            maximize(lambda_min(Sym(P*anyM*Q)))
            subject to
            norm(P)<=1;
        cvx_end
    end
