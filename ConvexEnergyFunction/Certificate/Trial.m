clc;
clearvars;
addpath('../QuadRelaxation/');
mpc     =   loadcase('case3sc.m');
[Ps,Qs] =   PQ(mpc);
res     =   ComputeSball('case3sc.m','case3scres.mat',.7,0*Ps,0*Qs);

