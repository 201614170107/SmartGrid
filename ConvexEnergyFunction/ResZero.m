function mpc=ResZero(mpc)
mpc.bus(:,5:6)=    0;
mpc.branch(:,4)=   (mpc.branch(:,4).^2+mpc.branch(:,3).^2)./(mpc.branch(:,4));
mpc.branch(:,3)=   0;
mpc.branch(:,5)=   0;
mpc.branch(:,9)=   0;
end