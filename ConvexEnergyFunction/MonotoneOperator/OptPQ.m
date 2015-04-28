function [vref,Popt,rhopt]=OptPQ(fname)
mpc             =   loadcase(fname);
define_constants;
fname           =   strcat(fname,'Popt.mat');
[Matsx,Matsy]   =   MakeMats(mpc);
n               =   size(Matsx,1);
m               =   size(Matsx,3);
sb              =   find(mpc.bus(:,2)==3);
vf              =   mpc.bus(sb,VM);
nsb             =   find(~(mpc.bus(:,2)==3));
th              =   linspace(0,2*pi,10);
thopt           =   0;
Popt            =   eye(n);
rhopt           =   0;
Ms              =   Matsx(:,:,nsb);
Ns              =   Matsy(:,:,nsb);

for it=1:length(th)
    V   =   vf*exp(1i*th(it));
    x   =   real(V);
    y   =   imag(V);
    M   =   sum(Matsx,3)*x+sum(Matsy,3)*y;
    [rhon,P]=OptP(M,Ms,Ns,rhopt);
            
    if(rhon>rhopt)
        rhopt   =   rhon;
        thopt   =   th(it);
        vref    =   vf*[cos(thopt);sin(thopt)];
        Popt    =   P;
        save(fname,'Popt','rhopt','vref');
    end
     disp(rhopt);
   
end




end

function [rhop,Pop]=OptP(M,Ms,Ns,rhopt)
    rhop    =   0;
    rhos    =   [rhopt,.5];
    n       =   size(M,1);
    Pop     =   eye(n);
   
    for kt=1:5        
        [chk,P] =   CheckRhoVlim(mean(rhos),M,Ms,Ns);
        if(chk)
            rhop    =   mean(rhos);
            rhos    =   [mean(rhos),max(rhos)];
            Pop     =   P;
        else
            rhos=[min(rhos),mean(rhos)];
        end
    end
end
