opt = spot_sdp_default_options();
opt.verbose = 1;
prog=spotsosprog;
x=msspoly('x',3);
prog = prog.withIndeterminate(x);
[prog,lam]=prog.newFree(1);
prog=prog.withSOS(lam-x'*x);
sol=prog.minimize(0,@spot_frlib,opt);

