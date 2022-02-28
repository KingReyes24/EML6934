function pdot = myode2(t,X0)

% prameters used to solve for theta
param.v    = X0(3);
param.lamx = X0(4);
param.lamy = X0(5);
param.lamv = X0(6);
param.g    = 10;

% initial guess for theta
theta_guess = 0;

% solve for theta 
options = optimoptions('fsolve','Display','off','TolFun',1e-8);
theta   = fsolve(@solveControl,theta_guess,options,param);
% integrate 
pdot = [
    param.v*sin(theta);
    param.v*cos(theta);
    param.g*cos(theta);
    0;
    0;
    -param.lamx*sin(theta)-param.lamy*cos(theta)
    ];

end