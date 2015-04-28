function P=FindP(mpc,Vs)
    [Matsx,Matsy]=MakeMats(mpc);
    k   =   size(Vs,2);
    n   =   size(Matsx,1);
    Ms  =   zeros(n,n,k);   
    for it=1:k
        Ms(:,:,it)  =   sum(bsxfun(@times,Matsx,permute(real(Vs(:,it)),[3,2,1])),3)...
                            +sum(bsxfun(@times,Matsy,permute(imag(Vs(:,it)),[3,2,1])),3);        
    end
    Sym=@(anyx) anyx+anyx';
    cvx_begin
    variable P(n,n);
    maximize(1)
    subject to
    for it=1:k
        Sym(P*Ms(:,:,it))-eye(n)==semidefinite(n);
    end
    cvx_end
end