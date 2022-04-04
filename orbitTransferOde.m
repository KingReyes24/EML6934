function pdot = orbitTransferOde(t,X,param)

mew = param.mew;
T   = param.T;
ve  = param.ve;

r         = X(1);
theta     = X(2);
vr        = X(3);
vtheta    = X(4);
m         = X(5);
lamr      = X(6);
lamtheta  = X(7);
lamvr     = X(8);
lamvtheta = X(9);
lamm      = X(10);

beta = atan2(lamvr,lamvtheta);

rdot         = vr;
thetadot     = vtheta/r;
vrdot        = vtheta^2/r - mew/r^2 + T*sin(beta)/m;
vthetadot    = -vtheta*vr/r + T*cos(beta)/m;
mdot         = -T/ve;
lamrdot      = lamtheta*vtheta/r^2 + lamvr*(vtheta^2/r^2 - 2*mew/r^3) - lamvtheta*vr*vtheta/r^2;
lamthetadot  = 0;
lamvrdot     = -lamr + lamvtheta*vtheta/r;
lamvthetadot = -lamtheta/r - 2*lamvr*vtheta/r + lamvtheta*vr/r;
lammdot      = lamvr*T*sin(beta)/m^2 + lamvtheta*T*cos(beta)/m^2;


pdot = [rdot;
        thetadot;
        vrdot;
        vthetadot;
        mdot;
        lamrdot;
        lamthetadot;
        lamvrdot;
        lamvthetadot;
        lammdot];
    
end