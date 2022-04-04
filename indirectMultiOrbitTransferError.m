function E = indirectMultiOrbitTransferError(z,param)
E       = [];
r0      = param.r0;
theta0  = param.theta0;
vr0     = param.vr0;
vtheta0 = param.vtheta0;
m0      = param.m0;

tf   = z(end);
zvec = [r0; theta0; vr0; vtheta0; m0; z(1:end-1)];
zvec = reshape(zvec,param.numStates,param.k);

time    = linspace(param.t0,tf,param.k+1);
options = odeset('reltol',1e-6);

for idx = 1:param.k
    X0    = zvec(:,idx);
    tspan = [time(idx) time(idx+1)];

    [t,p] = ode113(@orbitTransferOde,tspan,X0,options,param);

    s = length(t);
    
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
    
    
    if idx == param.k
        
        beta = atan2(lamvrf,lamvthetaf);

        Hf = lamrf*vrf + lamthetaf*vthetaf/rf + lamvrf*(vthetaf^2/rf - param.mew/rf^2 + param.T*sin(beta)/mf)...
            + lamvthetaf*(param.T*cos(beta)/mf - vthetaf*vrf/rf) + lammf*(-param.T/param.ve);

        E_temp = [ rf        - param.rf;
                   vrf       - param.vrf;
                   vthetaf   - param.vthetaf;
                   lamthetaf - param.lamthetaf;
                   lammf     - param.lammf;
                   Hf        - param.Hf ];

        E = [E; E_temp];  
    else
        pint   = reshape(p(end,:),[],1);
        E_temp = pint - zvec(:,idx+1);

        E = [E; E_temp]; 
    end
end
end