%% EML6934 Optimal Control
% Bonus 2 Problem 1
% This script solves the hyper sensitive optimal control problem
% by using the multiple shooting method
close all; clear all; clc
format longg
format compact

addpath(genpath('C:\Users\elias\Documents\UF_Classes\EML6934'))

% Known conditons
param.x0  = 1;
param.xf  = 1;
param.t0  = 0;

% Number of intervals
param.k = 2;

% List of different final times
param.tf = 10;

% create time intervals
t_inter = param.tf/param.k;
t_start = param.t0;
t_end   = param.t0+t_inter;
for idx = 1:param.k
    tvec{idx} = [t_start t_end];
    t_start = t_end;
    t_end   = t_end+t_inter;
end
param.time = tvec;

% Create conditions that need to be solved
lam0 = 0.4;
Pvec = 0.001*ones(1,2*(param.k-1));
P    = [lam0 Pvec];

%% solve optimality conditions 
fprintf('\nFinal Time = %f\n',param.tf);
options  = optimoptions('fsolve','Display','Iter','TolFun',1e-8);
f        = @(x)multiShoot(x,param);
solution = fsolve(f,P,options);

%% Integreate with solved conditions to get trajectory

d_tau = 2/param.k;
tau_start = -1;
for idx = 1:param.k
    tau_end = tau_start+d_tau;
    tau{idx} = [tau_start tau_end];
    tau_start = tau_end;
end

options = odeset('reltol',1e-8,'AbsTol',[1e-8 1e-8]);
for k_idx = 1:param.k
    if k_idx == 1       
        x   = param.x0;
        lam = solution(k_idx);
    else
        x   = solution(2*k_idx-2);
        lam = solution(2*k_idx-1);  
    end
    
%     [t,P] = ode45(@myode,[param.time{k_idx}(1) param.time{k_idx}(2)],[x lam],options);
    [t,P] = ode113(@myode,[tau{k_idx}(1) tau{k_idx}(2)],[x lam],options);
    p = 5*P;
    
    figure(1); hold on
    plot(t, P(:,1))
    xlabel('Time (s)')
    ylabel('X(t)')
    set(gcf,'color','white')
    set(gca,'fontweight','bold','fontsize',12)
    str1 = sprintf('Trajectory that Minimized the Objective Function (%d Intervals)',param.k);
    title(str1);
    grid minor
    
    figure(2); hold on;
    plot(t,-P(:,2))
    xlabel('Time (s)')
    ylabel('U(t)')
    str2 = sprintf('Control that Minimized the Objective Function (%d Intervals)',param.k);
    title(str2);
    set(gcf,'color','white')
    set(gca,'fontweight','bold','fontsize',12)
    grid minor
end
