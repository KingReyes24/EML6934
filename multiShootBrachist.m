function E = multiShootBrachist(P,param)

% create time intervals
t_inter = P(end)/param.k;
t_start = param.t0;
t_end   = param.t0+t_inter;
for idx = 1:param.k
    tvec{idx} = [t_start t_end];
    t_start = t_end;
    t_end   = t_end+t_inter;
end

theta = 0;
E = [];

for k_idx = 1:param.k
    if k_idx == 1
        x   = [param.x0 param.y0 param.v0];
        lam = P(1:3);
    else
        xStart   = 4+(k_idx-2)*param.numVars;
        lamStart = 7+(k_idx-2)*param.numVars;      
        x   = P(xStart:xStart+2);
        lam = P(lamStart:lamStart+2);
    end
    X0 = [x lam];
    options = odeset('reltol',1e-8,'AbsTol',[1e-8 1e-8 1e-8 1e-8 1e-8 1e-8]);

    [t,P_int]   = ode113(@myode2,[tvec{k_idx}(1) tvec{k_idx}(2)],X0,options);
    
    s = length(t);
    
    xf   = P_int(s,1);
    yf   = P_int(s,2);
    vf   = P_int(s,3);
    lamx = P_int(s,4);
    lamy = P_int(s,5);
    lamv = P_int(s,6);
    
    % solve for theta with final params
    theta_guess = theta;
    param2.lamx = lamx;
    param2.lamy = lamy;
    param2.lamv = lamv;
    param2.v    = vf;
    param2.g    = 10;
    
    options = optimoptions('fsolve','Display','off','TolFun',1e-8);
    [theta,fval,exitflag,output] = fsolve(@solveControl,theta_guess,options,param2);
    if exitflag ~= 1
        disp('not reliable');
    end
    H = lamx*vf*sin(theta) + (lamy*vf + lamv*param2.g)*cos(theta);
    
%     E = [ xf - param.xf
%         yf - param.yf
%         lamv
%         H + 1 ];
%     
    
    if k_idx == param.k
        E_temp = [xf - param.xf; 
                  yf - param.yf;
                  lamv;
                  H+1];
              
        E = [E; E_temp];
    else
        p_idx = 6*k_idx-2;
        E_temp = [xf - P(p_idx); 
                  yf - P(p_idx+1);
                  vf - P(p_idx+2); 
                  lamx - P(p_idx+3);
                  lamy - P(p_idx+4); 
                  lamv - P(p_idx+5)];
        E = [E; E_temp];
    end
    
end
end