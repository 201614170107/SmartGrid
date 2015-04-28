function [Popt,Qopt]=SetPQ(mpc,varargin)
[Matsx,Matsy]=MakeMats(mpc);
n=size(Matsx,1);
m=size(Matsx,3);
Mx=sum(Matsx,3);
My=sum(Matsy,3);
Pin=eye(n)/n;
Qin=eye(n)/n;
if(nargin>1)
    Pin=varargin{1};
    Qin=varargin{2};
end
th=linspace(0,2*pi,15);
topt=-inf;
Sym=@(anyx)(anyx+anyx');

for it=1:length(th)
    V   =   exp(1i*th(it));
    x   =   real(V);
    y   =   imag(V);
    M   =   Mx*x+My*y;
    P   =   Pin;
    Q   =   Pin;
    told=   -inf;
    for jt=1:10
        P   =   OptP(M,Q);
        Q   =   OptQ(M,P);
        tn  =   min(eig(Sym(P*M*Q)));
        if(tn-told<1e-5)
            break;
        else
            told=tn;
        end
    end
    
    t=tn;
    disp(t);
    if(topt<t)
        Popt=P;
        Qopt=Q;
        topt=t;
        thopt=th(it);
    end
end

    
vref=[cos(thopt);sin(thopt)];
end