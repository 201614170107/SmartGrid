clc;
clearvars;
define_constants;
mpc =   loadcase('case6ww');
[Mv,Mbr,JF,nn]    =   MakePolarJac2(mpc);
J                 =   JF(zeros(nn,1));
sol               =   SolvePolarSOS(Mv,Mbr,mpc,pi/10,1,.96);

%%
n=size(J,1);
cvx_begin
variable W(n,n);
maximize(lambda_min(W*J+(W*J)'))
subject to
sigma_max(W)<=1;
cvx_end