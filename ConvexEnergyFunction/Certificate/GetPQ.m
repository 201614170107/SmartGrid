function [P,Q]=GetPQ(mpc)
    res =   runopf(mpc);
    define_constants;
    V   =   res.bus(:,VM).*exp(1i*res.bus(:,VA)*pi/180);
    Y   =   makeYbus(mpc);
    S   =   V.*(Y*V);
    P   =   real(S);
    Q   =   imag(S);
end