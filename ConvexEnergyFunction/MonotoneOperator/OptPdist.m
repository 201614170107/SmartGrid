function [optv,P,Q]=OptPdist(Ms,Mo,lp)
n       =   size(Mo,1);
Sym     =   @(anyx) anyx+anyx';
SD      =   @(anyx) anyx-diag(diag(anyx));
k       =   size(Ms,3);

if(lp)
    A=[];
    Ms  =   cat(3,Ms,Mo);
    for it=1:size(Ms,3)
        Ai  =   zeros(n^2,n^2);
        sen =   [];
        Mi  =   Ms(:,:,it);
        for jt=1:n^2
            [ii,jj] =   ind2sub([n,n],jt);
            N       =   sparse(ii,1,1,n,1)*Mi(jj,:);
            Nd      =   N-diag(diag(N));
            Ai(:,jt)=   [diag(N)-sum(Nd,2);Nd(ones(n)-eye(n)>0)];
        end
        A=[A;sparse(Ai)];
    end
    model=struct();
    model.A=sparse(A);
    model.obj=zeros(n^2,1);
    model.sense='>';
    model.rhs=1e-5*ones(size(A,1),1);
    result=gurobi(model);
    optv=strcmp(result.status,'OPTIMAL');
    if(optv>0)
        P=reshape(result.x,[n,n]);
        for it=1:size(Ms,3)
            [~,chk]=chol(Sym(P*Ms(:,:,it)));
            
        end
    else
        P=[];
    end
else
    [U,S,V] =   svd(Mo);
    for it=1:k
        Ms(:,:,it)  =   U'*(Ms(:,:,it)+Mo)*V;
    end
    Mo      =   S;
    
    cvx_begin
    variable P(n,1);
    maximize(1)
    subject to
    for it=1:size(Ms,3)
        Sym(diag(P)*sparse(Ms(:,:,it)))-1e-4*eye(n)==semidefinite(n);
    end
    Sym(diag(P)*sparse(Mo))-1e-4*eye(n)==semidefinite(n);
    cvx_end
    optv=cvx_optval>.1;
    
    P=diag(P)*U';
    Q=V;
end
end