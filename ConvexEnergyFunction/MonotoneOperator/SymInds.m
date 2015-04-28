function S=SymInds(i,j,x,n)
    S=(sparse(i,j,x,n,n)+sparse(j,i,x,n,n));
   % S=S-diag(diag(S))/2;
end