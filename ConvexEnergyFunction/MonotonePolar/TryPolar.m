clc;
clearvars;
addpath('../QuadRelaxation');
define_constants;
mpc =   loadcase('case14');
[Mats,JF,nn]      =   MakePolarJac(mpc);
%J                   =   JF(zeros(nn,1));
opt               =   SolveMoment(Mats,mpc,.58,0);
%[Qse,Qbr,V]           =   FormQuads2(mpc,.7);
%sol                 =   SolvePolarSOS(Mv,Mbr,mpc,1,.517);

