function [chk,sols]=QuadEqns2(Msx,Qsx,Qsw,rho)

nn  =   size(Msx,1);
rs  =   @(anyx) reshape(anyx,[size(anyx,1)*size(anyx,2),size(anyx,3)]);
rss =   @(anyx,anyn) reshape(anyx,anyn,anyn);



%Independent Variables%
x       =   sdpvar(nn,1);
z       =   sdpvar(nn,1);
s       =   sdpvar(size(Qsx,3),1);
vars    =   [x;s;z];

%%Setup Main Constraint Quadratics %%

Qx  = [];
for it=1:size(Qsx,3)
    Qx  =   [Qx;([1;x]'*Qsx(:,:,it)*[1;x]-s(it))];
end

Qw  = [];
for it=1:size(Qsw,3)
    Qw  =   [Qw;([1;x]'*Qsw(:,:,it)*[1;x])];
end

Qs  =   rho-s'*s;

%%Setup Multipliers for non-zero determinant condition %%
T               =   [];
Fs              =   [];
[T,Fs,lamx]     =   AddPoly(T,Fs,vars,false,size(Qx,1),1);
[T,Fs,lamw]     =   AddPoly(T,Fs,vars,true,size(Qw,1),1);
[T,Fs,lams]     =   AddPoly(T,Fs,vars,true,size(Qs,1),1);
[T,Fs,lam]     =    AddPoly(T,Fs,vars,false,1,1);

%%Setup Non-zero determinant condition %%

TP=[];
pol=[];
if(MultSimple)
TP              =   sdpvar(nn,nn,size(Msx,3));
pol             =   z'*(rss(rs(Msx)*[1;x],nn)'*rss(rs(TP)*[1;x],nn))*z;
else
TP              =   sdpvar(1,1);
Fs              =   [Fs,sos(TP)];
pol             =   TP*z'*(rss(rs(Msx)*[1;x],nn)'*rss(rs(Msx)*[1;x],nn))*z;
end
Fs              =   [Fs,sos(pol-lamw'*Qw-lamx'*Qx-lams'*Qs-lam*(z'*z-1)-1)];

[sol,~,~,res]   =   solvesos(Fs,0,options,[T(:);TP(:)]);

%%Setup Quadratic Multipliers for each constraint and enforce necessary SOS %%
chk     =   false*ones(length(Qw)+1,1);
sols    =   cell(length(Qw)+1,1);

for it=1:size(Qsw,3)
    T               =   [];
    Fs              =   [];
    [T,Fs,lamx]     =   AddPoly(T,Fs,vars,false,size(Qx,1),1);
    [T,Fs,lamw]     =   AddPoly(T,Fs,vars,true,size(Qw,1),1);
    [T,Fs,lams]     =   AddPoly(T,Fs,vars,true,size(Qs,1),1);
    [T,Fs,lam]      =   AddPoly(T,Fs,vars,true,1,1);
  %  Fs              =   [Fs,sos(lam*Qw(it)-lamx'*Qx-lamw'*Qw-lams'*Qs-1)];
    [sol,~,~,res]   =   solvesos(Fs,0,options,T(:));%,[],options);
    chk(it)         =   all(res<1e-5);
    sols{it}        =   sol;
end



chk(end)        =   all(res<1e-5);
sols{end}       =   sol;
end

function [T,Fs,pols]=AddPoly(T,Fs,vars,SOS,n,m)
pols    =   sdpvar(n,m);
for it=1:n
    for jt=1:m
        [pol,cc]    =   polynomial(vars,2);
        T           =   [T;cc(:)];
        pols(it,jt) =   pol;
        if(SOS)
            Fs  =   [Fs,sos(pol)];
        end
    end
end
end