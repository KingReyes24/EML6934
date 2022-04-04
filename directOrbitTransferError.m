function [Eineq, E] = directOrbitTransferError(z,param)
Eineq = [];
t0      = param.t0;
tf      = z(end); 
c       = z(1:end-1); 
X0      = [param.r0; param.theta0; param.vr0; param.vtheta0; param.m0];
% tspan   = [t0 tf];
tspan   = [-1 1];
options = odeset('reltol',1e-6);

[t,p] = ode113(@directOrbitTransferOde,tspan,X0,options,c,param,tf);

s = length(t);

rf      = p(s,1);
vrf     = p(s,3);
vthetaf = p(s,4);


E= [ rf      - param.rf;
     vrf     - param.vrf;
     vthetaf - param.vthetaf ];

end