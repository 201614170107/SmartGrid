clc;clearvars;
define_constants;
%%
addpath('Cases/');
A       =   load('caselist.mat');
A       =   A.A;
bs      =   zeros(length(A),1);

for it=1:length(A)
    casename=strcat(A{it},'lossless');
    mpc=loadcase(casename);
res=runpf(mpc);
p=checkConv(res,mpc);

if(p==1)
    fprintf('\n Condition True \n');
    bs(it)=1;
else
    fprintf('\n Condition False \n');    
end
end

%%
mpc=loadcase('case118lossless');
mpc.branch(:,9)=0;
mpc.bus(:,6)=0;
res     =   runpf(mpc);
PQ      =   mpc.bus(:,BUS_TYPE)==1;
PV      =   PQ|(mpc.bus(:,BUS_TYPE)==2);
n       =   sum(PQ)+sum(PV);
rho     =   log(res.bus(:,VM));
theta   =   res.bus(:,VA)*pi/180;
%%
alphas=linspace(3.6,3.8,40);
mpopt=mpoption;
mpopt.sdp_pf.ndisplay=10000;
nosol=zeros(length(alphas),1);
normg=zeros(length(alphas),1);
for it=1:length(alphas)
   mpc2=ScaleLoads(mpc,alphas(it),1); 
   func=minEnergy(mpc2);
   thopt   =   minFunc(func,zeros(n,1),struct('display','off'));
   [~,g]=func(thopt);
   normg(it)=norm(g);
   nosol(it)=insolvablepf(mpc2,mpoption);
end
%%
plot(alphas,normg,'-k',alphas,nosol*max(normg),'--r');
ApplyLabel('\kappa','xlabel',16);
ApplyLabel('\delta=1 ','title',16);
legend('Normalized ||\nabla E||','PF Insolvability','FontSize',14);
legend('boxoff');

%%
func    =   minEnergy(res);

thopt   =   minFunc(func,zeros(n,1));

rhoe    =   rho;
rhoe(PQ)=   thopt(1:sum(PQ));
the     =   theta;
the(PV) =   thopt(sum(PQ)+1:end);
%%
alphas=linspace(1,2.43,30);
del=alphas;
for it=1:length(alphas)
    del(it)=maxVD(mpc,alphas(it));
end
plot(del,alphas);

