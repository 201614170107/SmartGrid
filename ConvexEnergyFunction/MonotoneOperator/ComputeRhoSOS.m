function [chk,sol]=ComputeRhoSOS(Mx,Qse,Qsw)
nn          =   size(Mx,1);
rs          =   @(anyx) reshape(anyx,[size(anyx,1)*size(anyx,2),size(anyx,3)]);
rss         =   @(anyx,anyn) reshape(anyx,anyn,anyn);
vec         =   @(anyx) anyx(:);
prog        =   spotsosprog;
x           =   msspoly('x',nn);
z           =   msspoly('z',nn);
prog        =   prog.withIndeterminate([x;z]);

%{
for it=2:size(Mx,3)
    Mx(:,:,it)=(Mx(:,:,1)\Mx(:,:,it));
end

Mx(:,:,1)   =   speye(nn);
%}

%Feron       =   2;

%if(Feron==2)
%    [prog,TP]   =   prog.newFree(nn*nn*size(Mx,3));
%    TP          =   reshape(TP,[nn*nn,size(Mx,3)]);    
%    pol         =   z'*(rss(TP*[1;x],nn)*rss(rs(Mx)*[1;x],nn)-eye(nn))*z;
%elseif(Feron==1)
%    [prog,TP]   =   prog.newFree(nn*nn);
%    TP          =   reshape(TP,[nn,nn]);
%    pol         =   z'*(TP*rss(rs(Mx)*[1;x],nn)-eye(nn))*z;
%else
%    [prog,TP]   =   prog.newFree(nn*nn);
%    TP          =   reshape(TP,[nn,nn]);
%    pol         =   z'*(rss(rs(Mx)*[1;x],nn)'*TP*rss(rs(Mx)*[1;x],nn)-eye(nn))*z;
%end
[prog,Me]   =   makePolyMat(prog,[x;z],size(Qse,3),1,0:2,-1);
[prog,Mw]   =   makePolyMat(prog,[x;z],size(Qsw,3),1,0:2,2);
[prog,Mm]   =   makePolyMat(prog,[x;z],nn,1,0:2,-1);



Qe  =   [];
for it=1:size(Qse,3)
    Qe=[Qe;[1;x]'*Qse(:,:,it)*[1;x]];
end

Qw  =   [];
for it=1:size(Qsw,3)
    Qw=[Qw;[1;x]'*Qsw(:,:,it)*[1;x]];
end


prog                =   prog.withSOS(Mm'*((rss(rs(Mx)*[1;x],nn))*z)-Me'*Qe-Mw'*Qw-z'*z);
options             =   spot_sdp_default_options();
options.verbose     =   1;
sol                 =   prog.minimize(0,@spot_frlib,options);
chk                 =   spotprogsol.statusIsPrimalFeasible(sol.status);


end

function [pr,laml]=makePolyMat(pr,vars,n,m,ds,SOSp)

if(SOSp==2)
    [pr,laml]   =   pr.newSOSPoly(monomials(vars,ds),n*m);
elseif(SOSp==1)
    [pr,laml]   =   pr.newSDSOSPoly(monomials(vars,ds),n*m);
elseif(SOSp==0)
    [pr,laml]   =   pr.newDSOSPoly(monomials(vars,ds),n*m);
else
    [pr,laml]   =   pr.newFreePoly(monomials(vars,ds),n*m);
end
laml        =   reshape(laml,[n,m]);

end
