function [P,Q,Plims,Qlims,V]=GetPQ(mpc)
    res     =   runpf(mpc);
    define_constants;
    Ps   =   -mpc.bus(:,PD);
    Qs   =   -mpc.bus(:,QD);
    Pmin =  Ps(mpc.gen(:,1))+mpc.gen(:,PMIN);
    Pmax =  Ps(mpc.gen(:,1))+mpc.gen(:,PMAX);
    Qmin =  Qs(mpc.gen(:,1))+mpc.gen(:,QMIN);
    Qmax =  Qs(mpc.gen(:,1))+mpc.gen(:,QMAX);
    Plims=  [Pmin,Pmax]/mpc.baseMVA;
    Qlims=  [Qmin,Qmax]/mpc.baseMVA;
    V   =   res.bus(:,VM).*exp(1i*res.bus(:,VA)*pi/180);
    Y   =   makeYbus(mpc);
    S   =   V.*conj(Y*V);
    P   =   real(S);
    Q   =   imag(S);
end