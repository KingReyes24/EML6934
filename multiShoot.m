function E = multiShoot(P,param)

options = odeset('reltol',1e-8,'AbsTol',[1e-8 1e-8]);

E = [];
% create time intervals
% d_tau = 2/param.k;
% tau_start = -1;
% for idx = 1:param.k
%     tau_end = tau_start+d_tau;
%     tau{idx} = [tau_start tau_end];
%     tau_start = tau_end;
% end

for k_idx = 1:param.k
    if k_idx == 1       
        x   = param.x0;
        lam = P(k_idx);
    else
        x   = P(2*k_idx-2);
        lam = P(2*k_idx-1);  
    end
% %     tspan = [-1 1];
  [t,P_int] = ode113(@myode,[param.time{k_idx}(1) param.time{k_idx}(2)],[x lam],options);    
%     [t,P_int] = ode113(@myode,[tau{k_idx}(1) tau{k_idx}(2)],[x lam],options);    
    
    size_t    = numel(t);
        
    if k_idx == param.k
        E_temp = P_int(size_t,1) - param.xf;
        E = [E; E_temp];
    else
        E_temp = [P_int(size_t,1) - P(2*k_idx); P_int(size_t,2) - P(2*k_idx+1)];
        E = [E; E_temp];
    end
end

end
