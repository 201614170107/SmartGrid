%{
function [chk,T]=QuadEqns(Msx,Qsx,Qsw,rho)
    Sym=@(anyx) (anyx+anyx');
    nn  =   size(Msx,1);
    n   =   size(Qsx,1);
    rs  =   @(anyx) reshape(anyx,[size(anyx,1)*size(anyx,2),size(anyx,3)]);
    sq  =   @(anyx,anyi,anyj) squeeze(anyx(anyi,anyj,:));
    rss =   @(anyx,anyn) reshape(anyx,[anyn,anyn,size(anyx,2)]);
    blZ =   @(anyx,anyi,anyj,sz) anyx((anyi-1)*sz+(1:sz),(anyj-1)*sz+(1:sz));
    
    cvx_begin
        variable P(nn,nn,size(Msx,3))
        variable T(nn,nn,size(Qsx,3));
        variable X(nn,nn,size(Qsx,3));
        variable U(nn,nn,size(Qsw,3));
        variable Z(nn*n,nn*n) symmetric;
        variable lambda(size(Qsx,3),size(Qsw,3));
        variable mul(size(Qsw,3),size(Qsw,3));
        variable kappa(size(Qsw,3));
        variable t(size(Qsw,3));
        maximize(1)
        subject to
        
        for ii=1:n
            for jj=(ii+1):n
                Sym(P(:,:,ii)*(squeeze(Msx(:,:,jj))))+...
                        Sym(P(:,:,jj)*(squeeze(Msx(:,:,ii))))-2*rss(rs(T)*sq(Qsx,ii,jj),nn)+...
                            -2*rss(rs(U)*sq(Qsw,ii,jj),nn)==2*blZ(Z,ii,jj,nn);
            end
        end
        for ii=2:n
           Sym(P(:,:,ii)*(squeeze(Msx(:,:,ii))))-rss(rs(T)*sq(Qsx,ii,ii),nn)...
                            -rss(rs(U)*sq(Qsw,ii,ii),nn)-blZ(Z,ii,ii,nn)==0;%-2*Z((ii-1)*nn+1:ii*nn,(jj-1)*nn+1:jj*nn)
        end
        
        for ii=1:size(T,3)
            T(:,:,ii)==T(:,:,ii)';
            X(:,:,ii)-rho*T(:,:,ii)==semidefinite(nn);
            X(:,:,ii)+rho*T(:,:,ii)==semidefinite(nn);
        end
        
        Sym(P(:,:,1)*(squeeze(Msx(:,:,1))))-rss(rs(T)*sq(Qsx,1,1),nn)-sum(X,3)...
               -rss(rs(U)*sq(Qsw,1,1),nn)...
                -blZ(Z,1,1,nn)-eye(nn)==semidefinite(nn);
        %T==0;
        for ii=1:size(U,3)
            U(:,:,ii)==semidefinite(nn);
        end
        
        for ii=1:size(Qsw,3)
            Qsw(:,:,ii)*kappa(ii)-rss(rs(Qsx)*lambda(:,ii),n)-rss(rs(Qsw)*mul(:,ii),n)-t(ii)*sparse(1,1,1,n,n)-eye(n)==semidefinite(n);
            t(ii)>=norm(lambda(:,ii),1)*rho;
        end
        mul>=0;
        kappa>=0;
        Z==semidefinite(n*nn);
        %Z(1:nn,(nn+1):n*nn)==0;
    cvx_end
    chk=(cvx_optval>0);
    T=cat(3,Sym(P(:,:,1)*(squeeze(Msx(:,:,1))))-rss(rs(T)*sq(Qsx,1,1),nn)...
               -rss(rs(U)*sq(Qsw,1,1),nn)...
                -blZ(Z,1,1,nn),T);
end
%}