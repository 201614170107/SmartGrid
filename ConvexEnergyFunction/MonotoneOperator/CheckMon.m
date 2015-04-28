function chk=CheckMon(mpc,P,sol)
    solx            =   real(sol);
    soly            =   imag(sol);
    [Matsx,Matsy]   =   MakeMats(mpc);
    M               =   sum(bsxfun(@times,Matsx,permute(solx,[3,2,1])),3)+...
                            sum(bsxfun(@times,Matsy,permute(soly,[3,2,1])),3);
    Sym             =   @(anyx) anyx+anyx';
    [~,p]           =   chol(Sym(P*M));
    chk             =   (p==0);
end