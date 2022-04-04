function [Eineq, E] = directMultiOrbitTransferError(z,param)
E       = [];
Eineq   = [];

% initial conditions
r0      = param.r0;
theta0  = param.theta0;
vr0     = param.vr0;
vtheta0 = param.vtheta0;
m0      = param.m0;
P0      = [r0; theta0; vr0; vtheta0; m0];

% seperate the unknowns 
% P is the unkown states
P_end  = param.numStates*(param.k-1);
P_tmp  = z(1:P_end);
P      = [P0;P_tmp];
P      = reshape(P,param.numStates,[]); % reshape to make each column an interval
% c is the uknown coeffients
c_end  = numel(z)-1;
c_list = z(P_end+1:c_end);
c_list = reshape(c_list,param.numCoeff,[]); % reshape to make each column an interval
tf     = z(end); 

% create grid to integrate on tau [-1 1]
tau     = linspace(-1,1,param.k+1);
options = odeset('reltol',1e-6);

for idx = 1:param.k
    % extract coefficents for the specific interval
    c  = c_list(:,idx);
    % extract states for the specific interval
    X0 = P(:,idx);
    % specific integartion span
    tspan = [tau(idx) tau(idx+1)];
    % solve for conditions at the end of the interval
    [t,p] = ode113(@directOrbitTransferOde,tspan,X0,options,c,param,tf);
    
    s = length(t);
 
    rf      = p(s,1);
    vrf     = p(s,3);
    vthetaf = p(s,4);
    
    % get the errors
    if idx == param.k
        % boundary conditions at final interval
        E_temp = [ rf        - param.rf;
                   vrf       - param.vrf;
                   vthetaf   - param.vthetaf ];

        E = [E; E_temp];  
    else
        % errors at each interval k-1
        pint   = reshape(p(end,:),[],1);
        E_temp = pint - P(:,idx+1);

        E = [E; E_temp]; 
    end
end
end