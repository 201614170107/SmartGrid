clc;
clearvars;
define_constants;
mpc =   loadcase('case14');
[Mv,Mbr,JF,nn]      =   MakePolarJac2(mpc);
J                   =   JF(zeros(nn,1));
opt                 =   SolveMoment(Mv,Mbr,mpc,.9,pi/4);
%sol                 =   SolvePolarSOS(Mv,Mbr,mpc,1,.517);

