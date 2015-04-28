clc;clearvars;
[vref,Popt,rhopt]=OptPQ('case3.m');

%%
mpc=loadcase('case3.m');
define_constants;
[Matx,Maty]=MakeMats(mpc);

th=atan(vref(2)/vref(1))+linspace(-pi/2.2,pi/2.2,5);
V=mpc.bus(2,VM);
[xs,ys]=ndgrid(vref(1)+linspace(-4,4,100),vref(2)+linspace(-4,4,100));

for it=1:length(th)
    chk=zeros(size(xs));
    for jt=1:numel(xs)
        Vx=[vref(1);V*cos(th(it));xs(jt)];
        Vy=[vref(2);V*sin(th(it));ys(jt)];
        M   =   sum(bsxfun(@times,Matx,permute(Vx,[3,2,1])),3)+...
                    sum(bsxfun(@times,Maty,permute(Vy,[3,2,1])),3);
        M   =   (Popt*M);
        [~,chk(jt)] =   chol(M+M');
        chk(jt)     =   chk(jt)==0;
    end
    x=xs(chk>0);x=x(:)-vref(1);
    y=ys(chk>0);y=y(:)-vref(2);
    K=convhull(x,y);
   figure(it);clf;  
    %plot(x,y,'.')
  
 axis equal
 hold on
 fill ( x(K), y(K), 'b','facealpha', 0.5 ); 
 %scatter(0,0,50,'filled');
 hold off
 ApplyLabel('$V^x$','xlabel',14);
 ApplyLabel('$V^y$','ylabel',14);
 z=(th(it)-atan(vref(2)/vref(1)))/pi;
 ApplyLabel(sprintf('$\\theta = %.2f \\pi$',z),'title',14);
  
end


%th=atan(vf(2)/vf(1));

%[xs,ys,zs]=ndgrid(linspace(vf(1)-.5,vf(1)+.5,50),linspace(vf(2)-.5,vf(2)+.5,50),linspace(th-pi/2,th+pi/2,50));

%[xs,ys]=ndgrid(linspace(vf(1)-2,vf(1)+2,50),linspace(vf(2)-2,vf(2)+2,50));

%vs=[xs(:),ys(:),zs(:)];
%chks=zeros(size(vs,1),1);
%%
define_constants;
res=runpf(mpc);
inj =   -res.bus(:,PD)-1i*res.bus(:,QD);
inj(res.gen(:,1))   =    inj(res.gen(:,1))+res.gen(:,PG)+res.gen(:,QG)*1i;

inj=inj/mpc.baseMVA;
V=res.bus(:,VM).*exp(1i*(res.bus(:,VA))*pi/180)*(vf(1)+vf(2)*1i);

%[fmin,CheckM]=minRegGap(mpc,vf,P,inj,res.bus(:,VM));
%%
xopt=minFunc(fmin,[vf(1)*ones(n-1,1);vf(2)*ones(n-1,1)]);
Ve=zeros(n,1);
Ve(sb)=vf(1)+1i*vf(2);
Ve(setdiff(1:n,sb))=xopt(1:n-1)+1i*xopt(n:2*(n-1));

Y=makeYbus(mpc);
S=inj-Ve.*conj(Y*Ve);



%%
for it=1:size(vs,1)
    x=[vf(1);cos(vs(it,3));vs(it,1)];
    y=[vf(2);sin(vs(it,3));vs(it,2)];
    M=sum(bsxfun(@times,Matx,permute(x,[3,2,1])),3)+...
            sum(bsxfun(@times,Maty,permute(y,[3,2,1])),3);
    M=M'*P+P*M;
    [~,p]=chol(M);
    chks(it)=(p==0);
  
end
%%
X=vs(chks>0,:);
X=[sqrt(X(:,1).^2+X(:,2).^2),(atan(X(:,2)./X(:,1))-th)*180/pi,(X(:,3)-th)*180/pi];

figure(1);clf;
%scatter3(X(:,1),X(:,2),X(:,3),20);
scatter(X(:,1),X(:,2),20);
ApplyLabel('V_3 (pu)','xlabel',14);
ApplyLabel('\theta_3 (degrees)','ylabel',14);

figure(2);clf;
scatter(X(:,2),X(:,3),20);
ApplyLabel('\theta_3 (degrees)','xlabel',14);
ApplyLabel('\theta_2 (degrees)','ylabel',14);


%%
