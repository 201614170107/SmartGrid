function [chk,sols]=QuadEqns(Msx,Qsx,Qsw,rho,varargin)

nn  =   size(Msx,1);
nm  =   size(Qsx,3);
rs  =   @(anyx) reshape(anyx,[size(anyx,1)*size(anyx,2),size(anyx,3)]);
rss =   @(anyx,anyn) reshape(anyx,anyn,anyn);

%   Solver Options     %
if(nargin<5)
    opts    =   QE_Default_Opt();
else
    opts    =   varargin{1};
end

if(opts.yalmip)
    solver          =   'frlib';
    options         =   sdpsettings('solver','frlib');
    opts.options    =   options;
else
    options         =   spot_sdp_default_options();
    options.verbose =   1;
    options.solver  =   @spot_frlib;
    opts.options    =   options;
end


% Do you want x-dependent multipliers a la Feron? %
%   Makes SOS program larger by considering x-depdendent multiplier %
MultFeron       =   opts.MultFeron;

mp              =   @(pr,xx,n,m) makePolyMat(pr,xx,n,m,2,0,opts);
mpS             =   @(pr,xx,n,m) makePolyMat(pr,xx,n,m,2,1,opts);


%Phase I (non-zero determinant) vs Phase II (invariance of non-zero determinant set)%
%Possibly less conserivative, but gives larger SOS problem%
mix             =   opts.mix;


%Independent Variables%
if(opts.yalmip)
    x   =   sdpvar(nn,1);
    z   =   sdpvar(nn,1);
    s   =   sdpvar(nm,1);
else
    x   = msspoly('x',nn);
    z   = msspoly('z',nn);
    s   = msspoly('s',nm);
end
%%Setup Main Constraint Quadratics %%

Qx  = [];
for it=1:nm
    Qx=[Qx;([1;x]'*Qsx(:,:,it)*[1;x]-s(it))];
end

Qw  = [];
for it=1:size(Qsw,3)
    Qw=[Qw;([1;x]'*Qsw(:,:,it)*[1;x])];
end

Qs  =   rho-s'*s;


%%Setup Multipliers for non-zero determinant condition %%
vars    =   [];
if(mix)
    vars        =   [x;z;s];
else
    vars        =   [x;z];
end
prog        =   getProg(opts);
prog        =   setIndeterminate(prog,vars,opts);
[prog,Mw]   =   mpS(prog,vars,size(Qsw,3),1);
[prog,Mm]   =   mp(prog,vars,1,1);


if(mix)
    [prog,Mx]   =   mp(prog,vars,nm,1);
    [prog,Ms]   =   mpS(prog,vars,size(Qs,1),1);
end

%%Setup Non-zero determinant condition %%

if(MultFeron)
    [prog,TP]   =   prog.newFree(nn*nn*size(Msx,3));
    TP          =   reshape(TP,nn*nn,size(Msx,3));
    pol         =   z'*(rss(rs(Msx)*[1;x],nn)'*rss(TP*[1;x],nn))*z;
else
    [prog,TP]   =   prog.newPSD(nn);
    pol         =   z'*(rss(rs(Msx)*[1;x],nn)'*TP*rss(rs(Msx)*[1;x],nn))*z;
end
if(mix)
    prog        =   AddSOS(prog,pol-Mw'*Qw-Mx'*Qx-Ms'*Qs-1,opts);
else
    prog        =   AddSOS(prog,pol-Mw'*Qw-Mm*(z'*z-1)-1,opts);
end

sols        =       cell(nm+1,1);
chk         =       zeros(nm+1,1);
[sols{1},chk(1)]      =       solveProg(prog,opts);

%%Setup invariance condition for constraints%%
vars        =   [x;s];

for it=1:size(Qsw,3)
    %If any one of the tests so far has failed, abort!%
    if(~all(chk(1:it)))
        break;
    end
    
    prog            =   spotsosprog;
    prog            =   prog.withIndeterminate(vars);
    
    %Parameters%
    [prog,lamx]     =    mp(prog,vars,size(Qsx,3),1);
    [prog,lamw]     =    mpS(prog,vars,size(Qsw,3),1);
    [prog,lam]      =    mpS(prog,vars,1,1);
    [prog,lams]     =    mpS(prog,vars,2,1);
    
    prog            =    AddSOS(prog,lam*Qw(it)-lamx'*Qx-lamw'*Qw-lams*Qs-1,opts);
    [sols{it+1},chk(it+1)]       =    solveProg(prog,opts);
end


end

function pr=getProg(opts)
if(opts.yalmip)
    pr      =   struct();
    pr.Fs   =   [];
    pr.Fv   =   [];
    pr.vars =   [];
else
    pr=spotsosprog;
end
end

function pr=setIndeterminate(pr,vars,opts)
if(opts.yalmip)
    pr.vars =   vars;
else
    pr  =   pr.withIndeterminate(vars);
end
end

function pr=AddSOS(pr,pol,opts)

if(opts.yalmip)
    Fs  =   pr.Fs;
    Fs  =   [Fs,sos(pol)];
    pr.Fs   =   Fs;
else
    if(opts.SOSpow==2)
        pr  =   pr.withSOS(pol);
    elseif(opts.SOSpow==1)
        pr  =   pr.withSDSOS(pol);
    else
        pr  =   pr.withDSOS(pol);
    end
end
end

function [sol,chk]=solveProg(prog,opts)
if(opts.yalmip)
    solvesos(prog.Fs,0,opts,options,prog.Fv(:));
else
    sol     =   prog.minimize(0,opts.options.solver,opts.options);
    chk     =   spotprogsol.statusIsPrimalFeasible(sol.status);
end
end


function [pr,laml]=makePolyMat(pr,vars,n,m,d,SOS,opts)
SOSp=opts.SOSpow;
if(opts.yalmip)
    pols    =   sdpvar(n,m);
    Fv      =   pr.Fv;
    Fs      =   pr.Fs;
    for ii=1:n
        for jj=1:m
            [poly,cc]    =   polynomial(vars,d);
            Fv          =   [Fv;cc(:)];
            pols(ii,jj) =   poly;
            if(SOS)
                Fs  =   [Fs,sos(pol)];
            end
        end
    end
    pr.Fv   =   Fv;
    pr.Fs   =   Fs;
    laml    =   pols;
else
    if(SOS)
        if(SOSp==2)
            [pr,laml]   =   pr.newSOSPoly(monomials(vars,0:d),n*m);
        elseif(SOSp==1)
            [pr,laml]   =   pr.newSDSOSPoly(monomials(vars,0:d),n*m);
        else
            [pr,laml]   =   pr.newDSOSPoly(monomials(vars,0:d),n*m);
        end
    else
        [pr,laml]   =   pr.newFreePoly(monomials(vars,0:d),n*m);
    end
    laml        =   reshape(laml,[n,m]);
end
end
