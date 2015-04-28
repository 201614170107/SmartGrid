clc;clearvars;
mpc=loadcase('case14');
[Matx,Maty]=MakeMats(mpc);
n=size(mpc.bus,1);
sb=find(mpc.bus(:,2)==3);
define_constants;
res=runpf(mpc);
inj =   -res.bus(:,PD)-1i*res.bus(:,QD);
inj(res.gen(:,1))   =    inj(res.gen(:,1))+res.gen(:,PG)+res.gen(:,QG)*1i;

inj=inj/mpc.baseMVA;
%[Popt,Qopt]=SetPQ(mpc);

[vf,Popt,Qopt]=findRef(mpc,0,mpc.bus(sb,VM));


[fmin,CheckM]=minRegGap(mpc,vf,Popt,Qopt,inj,res.bus(:,VM));

%%


V=res.bus(:,VM).*exp(1i*(res.bus(:,VA))*pi/180);%*(vf(1)+vf(2)*1i);
%V=V/V(sb);

n=(size(mpc.bus,1)-1);
xy=Qopt\[vf(1)*ones(n,1);vf(2)*ones(n,1)];
xyopt=minFunc(fmin,xy,struct('MaxIter',10000,'MaxFunEvals',10000));
xy=Qopt*xyopt;
Vres=zeros(n+1,1);
Vres(~(mpc.bus(:,2)==3))=(xy(1:n)+1i*xy(n+1:end));
Vres((mpc.bus(:,2)==3))=V(1);


