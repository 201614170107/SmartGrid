clc;
clearvars;
define_constants
mpc=loadcase(case3);
result=runpf(mpc);

P=result.bus(:,3);
Q=result.bus(:,4);
P(result.gen(:,1))=P(result.gen(:,1))-result.gen(:,2);
P=P/mpc.baseMVA;
Q=Q/mpc.baseMVA;
Bmat=imag(makeYbus(mpc));
Mat=Bmat-diag(diag(Bmat));
alph=sum(Mat(:,(mpc.bus(:,BUS_TYPE)==2)|(mpc.bus(:,BUS_TYPE)==3)),2)./sum(Mat,2);
delt=.8./(1-alph);


[i,j]   =   find(abs(Bmat)>0);


[th1,th2]   =   ndgrid(linspace(-pi/3,pi/3,100),linspace(-pi/3,pi/3,100));
thgrid      =   [th1(:),th2(:)]';
Vgrid       =   thgrid;
E           =   zeros(size(thgrid,2),1);
options     = optimoptions('fsolve','Jacobian','on','display','final');
for it=1:size(thgrid,2)
    disp(it);
    th  =   [0;thgrid(:,it)];
    Cmat=   full(Bmat.*sparse(i,j,cos(th(i)-th(j)),3,3));
    sf  =   solvefun(Cmat(2:3,:),Q(2:3),zeros(2,1));
    [Vs,~,ef]  =   fsolve(sf,zeros(2,1),options);
    Vs  =   [1;exp(Vs(:))];
    Vgrid(:,it)=  Vs(2:3);
    if(~(ef==1))
        E(it)   =   nan;
    else
        E(it)   =   P'*th+log(Vs(2:3))'*Q(2:3)-.5*Vs'*Cmat*Vs;
    end
end

E   =   reshape(E,size(th1));
HE  =   zeros(size(E));
del =   th1(2,1)-th1(1,1);

for i=3:size(th1,1)-2
    for j=3:size(th2,2)-2
        fprintf('\n %d %d \n',i,j);    
        H       =   zeros(2);
        H(1,1)  =   (16*(E(i+1,j)+E(i-1,j))-E(i-2,j)-E(i+2,j)-30*E(i,j))/12;   
        H(1,2)  =   (E(i+1,j+1)+E(i-1,j-1)-E(i+1,j-1)-E(i-1,j+1))/4;
        H(2,1)  =   H(1,2);
        H(2,2)  =   (16*(E(i,j+1)+E(i,j-1))-E(i,j-2)-E(i,j+2)-30*E(i,j))/12;                   
        HE(i,j) =   double(~any(isnan(H(:)))&(H(1,1)>=0)&(det(H)>=0));
    end
end

d=sum(Bmat,2);
M11=d(2)-Bmat(1,2)./Vgrid(1,:)./cos(th1(:)')-Bmat(2,3)*Vgrid(2,:)./Vgrid(1,:)./cos(th1(:)'-th2(:)');
M22=d(3)-Bmat(1,3)./Vgrid(2,:)./cos(th2(:)')-Bmat(2,3)*Vgrid(1,:)./Vgrid(2,:)./cos(th1(:)'-th2(:)');
M12=-Bmat(2,3)./cos(th1(:)'-th2(:)');

HET=    ((M11>0)&(M22>0)&(M11.*M22>M12.^2));

HET=reshape(HET,size(th1));
HET(1:2,:)=0;
HET(:,1:2)=0;
HET(end-1:end,:)=0;
HET(:,end-1:end)=0;
%%
figure(1);colormap jet;imagesc([60,-60],[-60,60],(2*reshape(HET,size(th1))+reshape(HE,size(th1)))');



