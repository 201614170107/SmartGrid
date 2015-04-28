function mpc2=ScaleLoads(mpc,fac,rat)
if(nargin<3)
    rat=1;
end
    define_constants;
    mpc2=mpc;
    mpc2.bus(:,PD)=fac*mpc.bus(:,PD);
    mpc2.bus(:,QD)=rat*fac*mpc.bus(:,QD);
    mpc2.gen(:,PG)=fac*mpc.gen(:,PG);
end