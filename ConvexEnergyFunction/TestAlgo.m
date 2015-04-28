clc;clearvars;
define_constants;
casename='case14new';
mpco =   runpf(ResZero(loadcase(casename)));

alpha=   0;
betas=  linspace(.8,.85,30);
errs =  zeros(3,length(betas));
tims =  zeros(3,length(betas));

for it=1:length(betas)
beta            =   betas(it);
Struct          =   Convertmpc(mpco,alpha,beta);
[fun,A,b,lb,ub] =   minEnergyFun(Struct);
Opt = opti('fun',@(x)fun(x,false),'grad',@(x)fun(x,true),'lin',...
            A,-full(b),full(b),'bounds',full(lb),full(ub),...
            'options',optiset('display','off'));
tic;        
[x,fval] = solve(Opt,zeros(length(Struct.P)+length(Struct.Q),1));
tims(3,it)=toc;

[Struct,mpc]    =   Convertmpc(mpco,0,beta);
tic;        
result  =   runpf(mpc,mpoption('verbose',0,'out.all',0));
tims(1,it)=toc;
tic;
result2 =   rundcpf(mpc,mpoption('verbose',0,'out.all',0));
tims(2,it)=toc;

x1  =   [(result.bus(:,VA)-result.bus(1,VA))*pi/180;log(result.bus(Struct.PQInds,VM))];
x1  =   min(max(x1,lb),ub);
x2  =   [(result2.bus(:,VA)-result2.bus(1,VA))*pi/180;log(result2.bus(Struct.PQInds,VM))];
x   =   [x(1:size(result.bus,1))-x(1);x(size(result.bus,1)+1:end)];
errs(:,it)  =   [EvalSol(Struct,x1);EvalSol(Struct,x2);EvalSol(Struct,x)];
disp(it);
end
%%
figure(1);
subplot(2,1,1);
plot(betas,errs(1,:),'-k',betas,errs(2,:),':b',betas,errs(3,:),'-.r')

ApplyLabel('Relative Error','ylabel',14);
ApplyLabel('Load Scaling Factor','xlabel',14);
%xlim(betas);


subplot(2,1,2);
plot(betas,tims(1,:),'-k',betas,tims(2,:),':b',betas,tims(3,:),'-.r')
%xlim([1,1.5]);

ApplyLabel('Solution time','ylabel',14);
ApplyLabel('Load Scaling Factor','xlabel',14);
