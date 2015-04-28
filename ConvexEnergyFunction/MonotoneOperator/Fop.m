function fun=Fop(mpc,v,vo)
    n=size(mpc.bus,1);
    Y=makeYbus(mpc);
    sb=find(mpc.bus(:,2)==3);
    pv=find(mpc.bus(:,2)==2);
    pq=find(mpc.bus(:,2)==1);
    [Matx,Maty]=MakeMats(mpc);
    
    function [f,grad]=funf(xy)
        x=zeros(n,1);
        y=x;
        x(setdiff(1:n,sb))=xy(1:(n-1));
        y(setdiff(1:n,sb))=xy(n:2*(n-1));
        x(sb)=vo(1);
        y(sb)=vo(2);
        V=x+1i*y;
        S=V.*conj(Y*V);
        
        F=[real(S);imag(S)];
        
        F(n+pv)=x(pv).^2+y(pv).^2;
        F([sb;n+sb])=[];
        
        f=v'*F;
        
        M=sum(bsxfun(@times,Matx,permute(x,[3,2,1])),3)+...
            sum(bsxfun(@times,Maty,permute(y,[3,2,1])),3);
        grad=M'*v; 
    end
fun=@funf;
end