function [chk,sol]=MomentQuadSDP(PolyObj,Qse,Qsw,QQ)
n           =   size(Qse{1},1)-1;

tic;
[is,js,spM] =   GetInds(Qse,Qsw);
[Yi,ny,kf]  =   ComputeIndices(is,js,n);
fprintf('\n Time Elapsed to construct model: %f \n',toc);

tic;
cvx_begin
variable y(ny,1);
minimize(0)
subject to
GetY        =   @(anyix,anyM) GetYconst(anyM,y,kf,anyix);

%   Constraints on first (n+1)x(n+1) block  of moment matrix %
y(Yi(1:n+1,1:n+1))==semidefinite(n+1);
y(kf([0;0;0;0]))==1;

%   QQ %
for it=1:length(QQ)
    trace(y(Yi(1:n+1,1:n+1))'*QQ{it})==0;
end

%   Higher order Constraints %

%   Equality Constraints %
for it=1:length(Qse)
    ix=FindClique(Qse{it},spM);
    GetY(ix,Qse{it})==0;
    trace(Qse{it}'*y(Yi(1:n+1,1:n+1)))==0;
end

%   Inequality Constraints %
for it=1:length(Qsw)
    ix=FindClique(Qsw{it},spM);
    GetY(ix,Qsw{it})==semidefinite(length(ix));
    trace(Qsw{it}'*y(Yi(1:n+1,1:n+1)))>=0;
end

%   PSD Constraints on Moment matrix%
for it=1:length(spM)
    [ix,jx] =   formPairs(spM{it});
    Mc      =   FormIndsClq(ix,jx,kf);
    if(~isempty(Mc))
        y(Mc)==semidefinite(size(Mc,1));
    end
end

fprintf('\n Time Elapsed to construct CVX problem: %f \n',toc);
cvx_end
chk=isinf(cvx_optval);
sol=y;
end


function [is,js,spM]=GetInds(Qse,Qsw)
n           =   size(Qse{1},1);
M           =   zeros(n);
for it=1:length(Qse)
    M   =   M+CF(Qse{it});
end
for it=1:length(Qsw)
    M   =   M+CF(Qsw{it});
end
M   =   spones(M+M');
spM =   ComputeDecomposition(M);

M   =   zeros(n);
for it=1:length(spM)
    M=M+sparse(spM{it},1,1,n,1)*sparse(spM{it},1,1,n,1)';
end
M                   =   double(M>0);
M(~triu(true(n)))   =   0;
[is,js]             =   find(M);
is                  =   is-1;
js                  =   js-1;
sqis                =   is==0;
is(sqis)            =   [];
js(sqis)            =   [];
is                  =   [zeros(n,1);is];
js                  =   [(0:(n-1))';js];
end

function M=CF(M)
ix=[find(sum(abs(M),2))];
M=spones(sparse(ix,1,1,size(M,1),size(M,2))*sparse(ix,1,1,size(M,1),size(M,2))');
end

function Mc=FindClique(M,Ms)
n=size(M,1);
ix=find(sum(abs(M),2));
lc=inf;
Mc=[];
ls=zeros(length(Ms),1);
for jt=1:length(Ms)
    ic=sparse(Ms{jt},1,1,n,1)>0;
    ls(jt)=sum(ic(ix));
    if(all(ic(ix)))
        if(lc>length(Ms{jt}))
            Mc=Ms{jt};
            lc=length(Mc);
        end
        return;
    end
end
[~,jt]=max(ls);
Mc=Ms{jt};
if(isempty(Mc))
    disp('trouble');
end
end

function Yc=GetYconst(Q,y,kf,ix)
FI      =   @(anyi,anyj) (FormInds(ix,anyi,anyj,kf));
[ic,jc] =   find(Q);
Yc      =   [];
in      =   (jc>=ic);
ic      =   ic(in);
jc      =   jc(in);
for kt=1:length(ic)
    i   =   ic(kt);
    j   =   jc(kt);
    res =   FI(i-1,j-1);
    if(isempty(res))
        Yc=[];
        return;
    end
    res =   y(res);
    if(i==j)
        if(isempty(Yc))
            Yc  =   Q(i,i)*res;
        else
            Yc  =   Yc+Q(i,i)*res;
        end
    else
        if(isempty(Yc))
            Yc  =   2*Q(i,j)*res;
        else
            Yc  =   Yc+2*Q(i,j)*res;
        end
    end
end
end

function M=FormInds(ix,i,j,kf)
M   =   zeros(length(ix));
for it=1:length(ix)
    for jt=it:length(ix)
        ind    =   kf([ix(it)-1;ix(jt)-1;i;j]);
        if(isempty(ind))
            M=[];
            return;
        end
        M(jt,it)    =   ind;
        M(it,jt)    =   ind;
    end
end
end

function M=FormIndsClq(ix,jx,kf)
M     =  zeros(length(ix));
for it=1:length(ix)
    for jt=it:length(ix)
        ind         =   kf([ix(it)-1;ix(jt)-1;jx(it)-1;jx(jt)-1]);
        if(isempty(ind))
            M=[];
            return;
        end
        M(jt,it)    =   ind;
        M(it,jt)    =   ind;
    end
end
end

function [ixx,jxx]=formPairs(ix)
M   =   ones(length(ix));
M(~triu(true(length(ix))))=0;
[ixx,jxx]=find(M);
ixx=ix(ixx);
jxx=ix(jxx);
end

%{
if(~usecvx)
    prog    =   spotprog;
    [prog,y]=   prog.newFree(ny);
    GetY    =   @(anyix,anyM) GetYconst(anyM,y,kf,anyix);
    
    
    prog    =   prog.withPSD(y(Yi(1:n+1,1:n+1)));
    prog    =   prog.withEqs(y(kf([0;0;0;0]))-1);
    
    prog    =   prog.withEqs(trace(y(Yi(1:n+1,1:n+1))'*QQ));
    
    
    for it=1:length(Qse)
        ix      =   FindClique(Qse{it},spM);
        prog    =   prog.withEqs(GetY(ix,Qse{it}));
        prog    =   prog.withEqs(trace(Qse{it}'*y(Yi(1:n+1,1:n+1))));
    end
    
    for it=1:length(Qsw)
        ix      =   FindClique(Qsw{it},spM);
        prog    =   prog.withPSD(GetY(ix,Qsw{it}));
        prog    =   prog.withPos(trace(Qsw{it}'*y(Yi(1:n+1,1:n+1))));
    end
    
    for it=1:length(spM)
        [ix,jx] =   formPairs(spM{it});
        Mc      =   FormIndsClq(ix,jx,kf);
        prog    =   prog.withPSD(y(Mc));
    end
    
    options             =   spot_sdp_default_options();
    options.verbose     =   1;
    sol                 =   prog.minimize(0,@spot_frlib,options);
    chk                 =   spotprogsol.statusIsPrimalFeasible(sol.status);
else
%}
