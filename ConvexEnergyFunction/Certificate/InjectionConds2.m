function chk=InjectionConds2(Qsx,Qsw,rho)
rs  =   @(anyx) reshape(anyx,[size(anyx,1)*size(anyx,2),size(anyx,3)]);
sq  =   @(anyx,anyi,anyj) squeeze(anyx(anyi,anyj,:));
rss =   @(anyx,anyn) reshape(anyx,[anyn,anyn,size(anyx,2)]);
blZ =   @(anyx,anyi,anyj,sz) anyx((anyi-1)*sz+(1:sz),(anyj-1)*sz+(1:sz));
n   =   size(Qsx,1);
chk =   zeros(size(Qsw,3),1);
x   =   sdpvar(n-1,1);
s   =   sdpvar(size(Qsx,3),1);
lam =   sdpvar(size(Qsx,3),1);
mu  =   sdpvar(size(Qsw,3)+1,1);
    
qsx =   [];
qsw =   [];
T   =   [];
zsx =   [];
zsw =   [];
for it=1:size(Qsx,3)
    qsx     =   [qsx;[1;x]'*Qsx(:,:,it)*[1;x]]-s(it);   
    [pp,cc] =   polynomial([x;s],2);
    zsx     =   [zsx;pp];
    T       =   [T;cc(:)];
end
for it=1:size(Qsw,3)
    qsw     =   [qsw;[1;x]'*Qsw(:,:,it)*[1;x]];   
    [pp,cc] =   polynomial([x;s],2);
    zsw     =   [zsw;pp];
    T       =   [T;cc(:)];
end
qsw     =   [qsw;rho-s'*s];   
[pp,cc] =   polynomial([x;s],2);
zsw     =   [zsw;pp];
T       =   [T;cc(:)];
[z,cc]  =   polynomial([x;s],2);
T       =   [T;cc(:)];
Fs      =   [];
for it=1:size(zsw,1)
    Fs  =   [Fs,sos(zsw(it))];
end
Fs      =   [Fs,sos(z)];


options = sdpsettings('sos.newton',1,'sos.congruence',1,'solver','frlib');
chk=[];
for i=1:size(Qsw,3)
    F=[Fs,sos(qsw(i)*z-zsx'*qsx-zsw'*qsw-1)];
    %t = sdpvar(1,1);
    %F=[sos(t)];
   % for it=1:size(Qsw,3)
    %    F=[F,sos(mu(it))];
   % end
   % F=[F,sos(qsw(i)*t-lam'*qsx-mu'*qsw-1)];
    [sol,v,Q,res] = solvesos(F,0,options,[T(:)]); 
    chk=[chk;all(abs(res)<1e-4)];
end
chk=all(chk);
end