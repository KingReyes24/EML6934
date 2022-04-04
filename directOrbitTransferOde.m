function pdot = directOrbitTransferOde(t,X,c,param,tf)

% control
beta = polyval(c,t);

% given constants
mew    = param.mew;
T      = param.T;
ve     = param.ve;

% states
r      = X(1);
theta  = X(2);
vr     = X(3);
vtheta = X(4);
m      = X(5);

% differential equations
rdot         = vr;
thetadot     = vtheta/r;
vrdot        = vtheta^2/r - mew/r^2 + T*sin(beta)/m;
vthetadot    = -vtheta*vr/r + T*cos(beta)/m;
mdot         = -T/ve;

pdot = [rdot;
        thetadot;
        vrdot;
        vthetadot;
        mdot];
    
% need to multiply by (tf-t0)/2 since I am on the tau [-1 1] scale
pdot = pdot * (tf - param.t0)/2;
end