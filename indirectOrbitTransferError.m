function E = indirectOrbitTransferError(z,param)
% error function for indirect shooting

lam0 = z(1:end-1);
tf   = z(end);

% initial conditons for ode. uses boundary conditions and unkowns co-states
X0   = [param.r0; param.theta0; param.vr0; param.vtheta0; param.m0; lam0];

options = odeset('reltol',1e-6);
tspan   = [param.t0 tf];

[t,p]   = ode113(@orbitTransferOde,tspan,X0,options,param);

s = length(t);

% propagated states at tf
rf         = p(s,1);
thetaf     = p(s,2);
vrf        = p(s,3);
vthetaf    = p(s,4);
mf         = p(s,5);
lamrf      = p(s,6);
lamthetaf  = p(s,7);
lamvrf     = p(s,8);
lamvthetaf = p(s,9);
lammf      = p(s,10);

% Control
beta = atan2(lamvrf,lamvthetaf);

% Hamiltonian
Hf = lamrf*vrf + lamthetaf*vthetaf/rf + lamvrf*(vthetaf^2/rf - param.mew/rf^2 + param.T*sin(beta)/mf)...
    + lamvthetaf*(param.T*cos(beta)/mf - vthetaf*vrf/rf) + lammf*(-param.T/param.ve);

% conditions
E = [ rf        - param.rf;
      vrf       - param.vrf;
      vthetaf   - param.vthetaf;
      lamthetaf - param.lamthetaf;
      lammf     - param.lammf;
      Hf        - param.Hf];

end