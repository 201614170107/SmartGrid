function chk=CheckMonN(mpc,Vs)
Y=makeYbus(mpc);
Ix=real(Y*Vs);
Iy=imag(Y*Vs);
sb=find(mpc.bus(:,2)==3);
mpc.bus(sb,:)=[];
pv=find(mpc.bus(:,2)==2);
pq=find(mpc.bus(:,2)==1);
npv=length(pv);
Y(sb,:)=[];
Y(:,sb)=[];
Vs(sb)=[];
Vx=real(Vs);
Vy=imag(Vs);
Ix(sb)=[];
Iy(sb)=[];
na=size(pq,1);
nb=size(pv,1);
n=na+nb;
Ya=Y(pq,pq);
Yb=Y(pq,pv);
Yc=Y(pv,pv);
%{
Mbig=[diag(Vx(pq)),diag(Vy(pq)),diag(Ix(pq)),diag(Iy(pq)),zeros(na,4*nb);...
    -diag(Vy(pq)),diag(Vx(pq)),diag(Iy(pq)),-diag(Ix(pq)),zeros(na,4*nb);...
    eye(2*na),[-imag(Ya),real(Ya);-real(Ya),-imag(Ya)],zeros(2*na,2*nb),[-real(Yb),imag(Yb);-imag(Yb),-real(Yb)];...
    [zeros(nb,4*na),diag(Vx(pv)),diag(Vy(pv)),diag(Ix(pv)),diag(Iy(pv))];...
    zeros(nb,2*na),imag(Yb'),real(Yb'),zeros(nb),eye(nb),imag(Yc),real(Yc);...
    10*[zeros(nb,4*na+2*nb),diag(Vx(pv)),diag(Vy(pv))];...
    zeros(nb,2*na),real(Yb'),-imag(Yb'),-eye(nb),zeros(nb),real(Yc),-imag(Yc)];
%}

Mbig=[diag(Vx(pq)),diag(Vy(pq)),diag(Ix(pq)),diag(Iy(pq)),zeros(na,2*nb);...
    -diag(Vy(pq)),diag(Vx(pq)),diag(Iy(pq)),-diag(Ix(pq)),zeros(na,2*nb);...
    eye(2*na),[-imag(Ya),real(Ya);-real(Ya),-imag(Ya)],[-real(Yb),imag(Yb);-imag(Yb),-real(Yb)];...
    10*[zeros(nb,4*na),diag(Vx(pv)),diag(Vy(pv))];...    
    [zeros(nb,2*na),diag(Vx(pv))*real(Yb')+diag(Vy(pv))*imag(Yb'),diag(Vy(pv))*real(Yb')-diag(Vx(pv))*imag(Yb'),diag(Ix(pv))+diag(Vx(pv))*real(Yc)+diag(Vy(pv))*imag(Yc),diag(Iy(pv))+diag(Vy(pv))*real(Yc)-diag(Vx(pv))*imag(Yc)]]; 
%Mbig([sb;n+sb],:)=[];
%Mbig(:,[2*n+sb;3*n+sb])=[];

Mbig=Mbig+Mbig';

[~,chk]=chol(Mbig);
chk=(chk==0);
end