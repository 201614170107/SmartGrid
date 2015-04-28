clf;clearvars;
b=1;
pt=0.25;
p=pt;
q=pt;

ts=linspace(-pi,pi,300);
vs=linspace(-5,.5,300);

[Xt,Xv]=meshgrid(ts,vs);
E=p*Xt+q*Xv-b*exp(Xv).*cos(Xt)+.5*b*exp(2*Xv);

fun=@(x)[p+b*exp(x(2))*sin(x(1));q-b*exp(x(2))*cos(x(1))+b*exp(2*x(2))];

funst=@(th)-2*p/b-sin(th).*cos(th)+sin(th).*sqrt(cos(th).^2-4*q/b);
funtt=@(th)-2*p/b-sin(th).*cos(th)-sin(th).*sqrt(cos(th).^2-4*q/b);

funs=@(th) (imag(funst(th))==0).*funst(th)+(abs(imag(funst(th)))>0)*1;
funt=@(th) (imag(funtt(th))==0).*funtt(th)+(abs(imag(funtt(th)))>0)*1;

vals=funs(ts);

vals(abs(imag(vals))>0)=inf;

valt=funt(ts);
valt(abs(imag(valt))>0)=inf;

figure(1);clf;
plot(ts,vals);
hold on;plot(ts,valt,'--r');
hold on;plot(ts,0*ts,':k');

sol=zeros(4,2);
sol(1,1)=fzero(funs,pi-.1);
sol(2,1)=fzero(funt,pi-.1);
sol(3,1)=fzero(funs,0);
sol(4,1)=fzero(funt,0);
sol=sol(sin(sol(:,1))<0,:);
sol(:,2)=log(-p./(b*sin(sol(:,1))));

figure(2);clf;
imagesc(ts,vs,E);colormap hot;
hold on;
scatter(sol(:,1),sol(:,2),100,[1,0,0;.1,1,0],'fill');

hold on;
plot(ts(abs(ts)<pi/2),log(1./(2*cos(ts(abs(ts)<pi/2)))),'--g','linewidth',20);


