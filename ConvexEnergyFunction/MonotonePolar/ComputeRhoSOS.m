function [chk,sol]=ComputeRhoSOS(Mx,Qse,Qsw)

nn          =   size(Mx,1);
m           =   size(Mx,3);
 
Qse         =   Cleanup(Qse);
Qsw         =   Cleanup(Qsw);
M           =   Cleanup(reshape(Mx,[nn^2,m]));
Mx          =   M{1};
rs          =   @(anyx) reshape(anyx,[size(anyx,1)*size(anyx,2),size(anyx,3)]);
rss         =   @(anyx,anyn) reshape(anyx,anyn,anyn);
vec         =   @(anyx) anyx(:);
Sym         =   @(anyx) anyx+anyx';
prog        =   spotsosprog;
x           =   msspoly('x',m-1);
z           =   msspoly('z',nn);
prog        =   prog.withIndeterminate([x;z]);


[prog,Me]   =   makePolyMat(prog,[z],length(Qse),1,2,-1);
[prog,Mw]   =   makePolyMat(prog,[z],length(Qsw),1,2,2);
[prog,W]    =   prog.newFree(1);
W           =   W*eye(nn);



Qe  =   [];
for it=1:length(Qse)
    Qe=[Qe;[1;x]'*Qse{it}*[1;x]];
end

Qw  =   [];
for it=1:length(Qsw)
    Qw=[Qw;[1;x]'*Qsw{it}*[1;x]];
end


prog                =   prog.withSOS(z'*Sym(W*(rss(Mx*[1;x],nn)))*z-Me'*Qe-Mw'*Qw-z'*z);
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

function Qsn=Cleanup(Qs)
    Qsn =cell(size(Qs,3),1);
    for it=1:size(Qs,3)
        M       =   Qs(:,:,it);
        M(abs(M)<max(abs(M(:)))*1e-5)  =   0;
        M       =   sparse(M);
        Qsn{it} =   M;
    end
end
