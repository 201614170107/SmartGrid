clc;
clearvars;
mpc=loadcase('case9.m');

[Matsx,Qse,Qsi,nn,TQ,ff]=MakeMats2(mpc,.25,1);
%%xy=minFunc(ff,randn(nn,1),struct('derivativeCheck','on'));

define_constants;
nsb=find(~(mpc.bus(:,BUS_TYPE)==3));
nrb=size(mpc.branch,1);
nb=size(mpc.bus,1);

[cp,sol]=ComputeRhoSOS(Matsx,Qse,Qsi);
%{
%%
load('case14lossless.mPopt.mat');
mpc=loadcase('case14lossless.m');
%fmin=minRegGap(mpc,vref,Popt);
n=size(mpc.bus,1)-1;
%xy=minFunc(fmin,[vref(1)*ones(n,1);vref(2)*ones(n,1)]);


define_constants;
mpopt=mpoption('verbose',0);

[fac,rat]   =   ndgrid(linspace(1,7,10),linspace(pi/6,pi/3,10));
resa        =   zeros(numel(fac),3);


for it=1:length(fac(:))

mpc2=ScaleLoads(mpc,fac(it),rat(it));
fmin=minRegGap(mpc2,vref,Popt);
%[xy,ff]=minFunc(fmin,[vref(1)*ones(n,1);vref(2)*ones(n,1)]);
res=runpf(mpc2,mpopt);
if(res.success)
Vs=res.bus(:,VM).*exp(1i*res.bus(:,VA)*pi/180)*(vref(1)+1i*vref(2))/norm(vref);
resa(it,1)=CheckMon(mpc2,Popt,Vs);
resa(it,2)=1;
resa(it,3)=CheckEF(mpc2,Vs);
else
resa(it,1)=0;    
resa(it,2)=1-insolvablepf(mpc2,mpopt);
resa(it,3)=0;
end
end

%Vss=[Vss(:,qs>0),Vss(sb)*ones(size(mpc.bus,1),1)];
%%
mpc=loadcase('case9.m');
[Mx,My]=MakeMats
%}
