%% EML6934 Optimal Control
%  Name:       Elias Reyes
%  Date:       04 April 2022
%  Assignment: Midterm
%  Goal:       Solve the orbit transfer problem using indirect shooting
close all; clear all; clc
format longg
format compact

% given constants 
param.T         = 0.1405;
param.mew       = 1;
param.ve        = 1.8658344;
% boundary conditions
param.r0        = 1;
param.theta0    = 0;
param.vr0       = 0;
param.vtheta0   = sqrt(param.mew/param.r0);
param.rf        = 1.5;
param.vrf       = 0;
param.vthetaf   = sqrt(param.mew/param.rf);
param.t0        = 0;
param.m0        = 1;
% transversality coniditions 
param.lamthetaf = 0;
param.lammf     = 1;
param.Hf        = 0; 
% initial guesses for unkown parameters 
lamr      = -2;
lamtheta  = 0;
lamvtheta = 2;
lamvr     = 2;
lamm      = 2;
tf        = 3.5;
zguess    = [lamr; lamtheta; lamvtheta; lamvr; lamm; tf];

options = optimoptions('fsolve','Display','Iter','TolFun',1e-8);

% solve the optimal control problem
f = @(x)indirectOrbitTransferError(x,param);
tic
solution =fsolve(f, zguess ,options);
toc

%%  Take optimal initial conditions and integrate them for plotting purposes
lam0 = solution(1:end-1); % optimal lambda values as t0
tf   = solution(end);     % optimal final time 
X0   = [param.r0; param.theta0; param.vr0; param.vtheta0; param.m0; lam0];

options = odeset('reltol',1e-6);
tspan   = [param.t0 tf];

% integrate with optimal conditions
[t,p]   = ode113(@orbitTransferOde,tspan,X0,options,param);

lamvrf     = p(:,8);
lamvthetaf = p(:,9);

% optimal control, use "mod" to ensure 0 -> 2pi
beta = mod(atan2(lamvrf,lamvthetaf),2*pi);

% plot figures
figure; hold on; grid minor
plot(t,p(:,1));
plot(t,p(:,2));
plot(t,p(:,3));
plot(t,p(:,4));
plot(t,p(:,5));
xlabel('t','Interpreter','LaTeX')
legend('r(t)','$\theta(t)$','vr(t)','v$\theta(t)$','m(t)','Interpreter','LaTeX')
set(gcf,'color','white')
set(gca,'fontweight','bold','fontsize',10)
title('States for trajectory that minimized fuel consumption')

figure; hold on; grid minor;
plot(t,beta*180/pi)
set(gcf,'color','white')
set(gca,'fontweight','bold','fontsize',10)
xlabel('t','Interpreter','LaTeX')
ylabel('$\beta(t)$','Interpreter','LaTeX')
title('Optimal Control that minimizes fuel consumption')
