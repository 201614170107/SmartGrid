clc;
clearvars;
mpc=loadcase('case3sc.m');

[Matsx,Matsw,Qssx,Qssw,fun,ff,TestQ,nn]=MakeMats(mpc,-.15,1);

define_constants;
nsb=find(~(mpc.bus(:,BUS_TYPE)==3));
nrb=size(mpc.branch,1);
nb=size(mpc.bus,1);

Matsx=Matsx(1:2*nb-2,1:2*nb-2,:);
[cp,sol]=QuadEqns(Matsx,Qssx,Qssw,.5);

