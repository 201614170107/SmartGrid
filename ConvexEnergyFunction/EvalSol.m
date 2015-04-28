function errnorm=EvalSol(Struct,x)
    n       =   size(Struct.P,1);
    m       =   size(Struct.Q,1);
    theta   =   x(1:n);
    V       =   zeros(n,1);
    V(Struct.PVInds)    =   Struct.Vs;
    V(Struct.PQInds)    =   exp(x(n+1:end));
    Bmat    =   Struct.Bmat;
    P       =   Struct.P;
    Q       =   Struct.Q;
    PQInds  =   Struct.PQInds;
    [i,j]   =   find(Bmat);
    Mat     =   Bmat.*sparse(i,j,cos(theta(i)-theta(j)),n,n).*(V*V');
    errnorm =   norm([P+sum(Bmat.*sparse(i,j,sin(theta(i)-theta(j)),n,n).*(V*V'),2);...
                      Q-sum(Mat(PQInds,:),2)])/norm([P;Q]);  
end