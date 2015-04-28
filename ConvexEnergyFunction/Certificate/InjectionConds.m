function chk=InjectionConds(Qsx,Qsw,rho)
rs  =   @(anyx) reshape(anyx,[size(anyx,1)*size(anyx,2),size(anyx,3)]);
sq  =   @(anyx,anyi,anyj) squeeze(anyx(anyi,anyj,:));
rss =   @(anyx,anyn) reshape(anyx,[anyn,anyn,size(anyx,2)]);
blZ =   @(anyx,anyi,anyj,sz) anyx((anyi-1)*sz+(1:sz),(anyj-1)*sz+(1:sz));
n   =   size(Qsx,1);
chk =   zeros(size(Qsw,3),1);
for i=1:size(Qsw,3)
    ch  =   [false;false];
    cvx_begin quiet
    variable lambda(size(Qsx,3),1);
    variable mul(size(Qsw,3),1);
    variable kappa(1);
    variable t(1)
    maximize(1)
    subject to
    rss(rs(Qsx)*lambda,n)-rss(rs(Qsw)*mul,n)+kappa*Qsw(:,:,i)-t*sparse(1,1,1,n,n)-eye(n)==semidefinite(n)
    mul>=0;    
    t>=rho*norm(lambda,1);
    cvx_end
    ch(1)   =   cvx_optval>0;
    
    cvx_begin quiet
    variable lambda(size(Qsx,3),1);
    variable mul(size(Qsw,3),1);
    variable kappa(1);
    variable t(1)
    maximize(1)
    subject to
    -(rss(rs(Qsx)*lambda,n)+rss(rs(Qsw)*mul,n)+kappa*Qsw(:,:,i)+t*sparse(1,1,1,n,n))-eye(n)==semidefinite(n)
    mul>=0;
    t>=rho*norm(lambda,1);
    cvx_end
    ch(2)   =   cvx_optval>0;
    
    chk(i)  =   ch(1)||ch(2);    
end
chk=all(chk);
end