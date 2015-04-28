function [fmin,Check]=minRegGap(mpc,vo,P)
[Mx,My] =   MakeMats(mpc);
inj     =   -mpc.bus(:,3)-1i*mpc.bus(:,4);
inj(mpc.gen(:,1))   =    inj(mpc.gen(:,1))+mpc.gen(:,2)+mpc.gen(:,3)*1i;
Sinj    =   inj/mpc.baseMVA;
vset    =   mpc.bus(:,8);
n       =   size(mpc.bus,1);
sb      =   find(mpc.bus(:,2)==3);
nsb     =   setdiff(1:n,sb);
Y       =   makeYbus(mpc);
pv      =   find(mpc.bus(:,2)==2);
Mscl    =   100*eye(2*(n-1));
R       =   chol(Mscl);
[Mo,Mall]  =    GetJacs();


    function [J,Jacs]=GetJacs()
        Jacs=zeros((2*(n-1))^2,2*(n-1));
        [~,J]=CheckM(zeros(2*(n-1)));
        J=J+J';
        for it=1:2*(n-1)
            [~,Jc]=CheckM(sparse(it,1,1,2*(n-1),1));
            Jc  =   Jc+Jc'-J;           
            Jacs(:,it)  =   Jc(:);
        end
    end

    function [res,gF]=CheckM(xy)
        [x,y]=makexy(xy);
        gF=sum(bsxfun(@times,Mx,permute(x,[3,2,1])),3)...
            +sum(bsxfun(@times,My,permute(y,[3,2,1])),3);
        gF=P*gF;
        [~,p]=chol(gF'+gF);
        res=(p==0);
    end

    function [x,y]=makexy(xy)
        x=zeros(n,1);
        y=x;
        x(nsb)=xy(1:(n-1));
        y(nsb)=xy(n:2*(n-1));
        x(sb)=vo(1);
        y(sb)=vo(2);
    end

    function F=formf(x,y)
        V   =   x+1i*y;
        S   =   V.*conj(Y*V);
        
        F   =   [real(S);imag(S)]-[real(Sinj);imag(Sinj)];
        
        F(n+pv) =   x(pv).^2+y(pv).^2-vset(pv).^2;
        F([sb;n+sb])    =   [];
        F   =   P*F;
    end

    function [f,grad]=funf(xy)
        [x,y]=makexy(xy);
        [res,gF]=CheckM(xy);
        if(~res)
            f=inf;
            grad=ones(size(xy));
        else
            F=formf(x,y);
            yp  =   Mscl\F;
            if(~CheckM(yp-xy))
                yp=ProjC(xy,yp,Mo,Mall,R);
            end
            f       =   F'*yp-.5*yp'*Mscl*yp;
            grad    =   F+(gF'-Mscl)*yp;
        end
        
    end




fmin=@funf;
Check=@CheckM;
end

function xyp=ProjC(xy,yp,M,Mall,R)
Sym     =   @(anyx) (anyx+anyx')/2;
n=size(xy,1);
cvx_begin quiet
variable xyp(n,1);
minimize(norm(R*(xyp-yp)))
subject to
M+Sym(reshape(Mall*(xy+xyp),size(M)))==semidefinite(size(M,1));
cvx_end
end