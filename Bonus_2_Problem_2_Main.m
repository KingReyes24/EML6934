%% EML6934 Optimal Control
% Bonus 2 Problem 2
% Solve the Brachistochrone problem using the multiple shooting method
close all; clear all; clc
format longg
format compact
addpath(genpath('C:\Users\elias\Documents\UF_Classes\EML6934'))
%%
% LOOK INTO k_idx =2. Always fails theta after awhile
%%
%known conditions
param.x0 = 0;
param.y0 = 0;
param.v0 = 0;
param.xf = 2;
param.yf = 2;
param.g  = 10;
param.t0 = 0;

% Number of intervals
param.k = 3;

% initial guesses for unkown parameters 
guess.lamx = 0;
guess.lamy = 0;
guess.lamv = 0;
guess.tf   = 1;

% Create conditions that need to be solved
param.numVars = 6;
lam0 = [-0.1 -0.1 -0.1];
Pvec = 0.001*ones(1,param.numVars*(param.k-1));
P    = [lam0 Pvec guess.tf];

options = optimoptions('fsolve','Display','Iter','TolFun',1e-8);

f = @(x)multiShootBrachist(x,param);
solution =fsolve(f,P,options);

%integrate with solution to create trajectory
t_inter = solution(end)/param.k;
t_start = param.t0;
t_end   = param.t0+t_inter;
for t_idx = 1:param.k
    tvec{t_idx} = [t_start t_end];
    t_start = t_end;
    t_end   = t_end+t_inter;
end

for k_idx = 1:param.k
    if k_idx == 1
        x   = [param.x0 param.y0 param.v0];
        lam = solution(1:3);
    else
        xStart   = 4+(k_idx-2)*param.numVars;
        lamStart = 7+(k_idx-2)*param.numVars;        
        x   = solution(xStart:xStart+2);
        lam = solution(lamStart:lamStart+2);
    end
    X0 = [x lam];
    options = odeset('reltol',1e-8,'AbsTol',[1e-8 1e-8 1e-8 1e-8 1e-8 1e-8]);
    
    [t,P_int]   = ode113(@myode2,[tvec{k_idx}(1) tvec{k_idx}(2)],X0,options);
    
    % solve for theta
    theta_guess = 0;
    theta = zeros(numel(t),1);
    for idx = 1:numel(t)
        param2.v    = P_int(idx,3);
        param2.lamx = P_int(idx,4);
        param2.lamy = P_int(idx,5);
        param2.lamv = P_int(idx,6);
        
        options = optimoptions('fsolve','Display','off','TolFun',1e-8);
        theta(idx) = fsolve(@solveControl,theta_guess,options,param2);
        theta_guess = theta(idx);
    end
    figure(1); hold on;
    plot(t,P_int(:,1),'b')
    plot(t,P_int(:,2),'r')
    plot(t,P_int(:,3),'g')
    plot(t([1 end]),P_int([1 end],1),'*k')
    plot(t([1 end]),P_int([1 end],2),'*k')
    plot(t([1 end]),P_int([1 end],3),'*k')
    str1 = sprintf('Trajectory that Minimized the Objective Function (%d Intervals)',param.k);
    title(str1);
    xlabel('Time (s)')
    legend('X(t)','Y(t)','V(t)','Interval Points')
    set(gcf,'color','white')
    set(gca,'fontweight','bold','fontsize',12)
    grid minor
    
    figure(2); hold on;
    plot(t,theta)
    xlabel('Time (s)')
    ylabel('Theta (deg)')
    str2 = sprintf('Control that Minimized the Objective Function (%d Intervals)',param.k);
    title(str2);
    set(gcf,'color','white')
    set(gca,'fontweight','bold','fontsize',12)
    grid minor
    
end