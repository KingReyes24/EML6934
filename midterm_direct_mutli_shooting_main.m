%% EML6934 Optimal Control
%  Name:       Elias Reyes
%  Date:       04 April 2022
%  Assignment: Midterm
%  Goal:       Solve the orbit transfer problem using direct multiple shooting
close all; clear all; clc
format longg
format compact

% specify number of intervals and desired degree polynomial degree
param.k         = 16;        % numbe of intervals
param.n         = 2;         % degree of polynomial used to parameterize the control

param.numCoeff  = param.n+1; % always 1 more coefficent then degree
param.numStates = 5;         % r,theta,vr,vtheta,m
%known conditions
param.T         = 0.1405;
param.mew       = 1;
param.ve        = 1.8658344;
% boundary conditions
param.r0        = 1;
param.theta0    = 0;
param.vr0       = 0;
param.vtheta0   = sqrt(param.mew/param.r0);
param.m0        = 1;
param.t0        = 0;
param.rf        = 1.5;
param.vrf       = 0;
param.vthetaf   = sqrt(param.mew/param.rf);
% initial guesses for unkown parameters 
tf       = 10;
pguess   = ones(param.numStates*(param.k-1),1);
cguess   = zeros(param.numCoeff*param.k,1);
tfguess  = tf;
zguess   = [pguess; cguess; tfguess];

% set up bounds
pmin  = -5*ones(numel(pguess),1);
pmax  =  5*ones(numel(pguess),1);
cmin  = -50*ones(param.numCoeff*param.k,1);
cmax  =  50*ones(param.numCoeff*param.k,1);
tfmin = 0;
tfmax = 11;
zmin  = [pmin; cmin; tfmin];
zmax  = [pmax; cmax; tfmax];
A     = [];
B     = [];
Aeq   = [];
Beq   = [];

options = optimset('MaxFunEvals',100000,'MaxIter',1000,'Display','Iter','TolFun',1e-4);

tic
solution = fmincon(@orbitTransferObj,zguess,A,B,Aeq,Beq,zmin,zmax,@directMultiOrbitTransferError,options,param);
toc

%% plot results
color   = ['b','r','g','c','m'];
r0      = param.r0;
theta0  = param.theta0;
vr0     = param.vr0;
vtheta0 = param.vtheta0;
m0      = param.m0;
P0      = [r0; theta0; vr0; vtheta0; m0];

% seperate the unknowns 
P_end  = param.numStates*(param.k-1);
P_tmp  = solution(1:P_end);
P      = [P0;P_tmp];
P      = reshape(P,param.numStates,[]);
c_end  = numel(solution)-1;
c_list = solution(P_end+1:c_end);
c_list = reshape(c_list,param.numCoeff,[]);
tf     = solution(end); 
t0     = param.t0;
% number of coefficients for polynomial

tau     = linspace(-1,1,param.k+1);
options = odeset('reltol',1e-6);

for idx = 1:param.k
    % extract coefficents for the specific interval
    c  = c_list(:,idx); 
    X0 = P(:,idx);

    tspan = [tau(idx) tau(idx+1)];
    [t,p] = ode113(@directOrbitTransferOde,tspan,X0,options,c,param,tf);
    
    % optimal control
    beta =  mod(polyval(c,t),2*pi);
    
    figure(1); hold on;
    for p_idx = 1:param.numStates
        plot(t,p(:,p_idx),color(p_idx));
        h =  plot(t([1 end]),p([1 end],p_idx),'*k');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    xlabel('$\tau$','Interpreter','LaTeX')
    legend('$r(\tau)$','$\theta(\tau)$','$v_r(\tau)$','$v_\theta(\tau)$','$m(\tau)$','Interpreter','LaTeX')
    set(gcf,'color','white')
    set(gca,'fontweight','bold','fontsize',10,'XMinorGrid','on','YMinorGrid','on')
    str1 = sprintf('States for trajectory that minimized fuel consumption (%d Intervals)',param.k);
    title(str1);
    
    figure(2); hold on; 
    plot(t([1 end]),beta([1 end])*180/pi,'*k')
    plot(t,beta*180/pi)
    set(gcf,'color','white')
    set(gca,'fontweight','bold','fontsize',10,'XMinorGrid','on','YMinorGrid','on')
    xlabel('$\tau$','Interpreter','LaTeX')
    ylabel('$\beta(\tau)$','Interpreter','LaTeX')
    str1 = sprintf('Optimal Control that minimizes fuel consumption (%d Intervals)',param.k);
    title(str1);
end
tf 
p(end,5)
c_list
% figure(1)
% print -depsc directStatesK16Poly5.eps;
% figure(2)
% print -depsc directControlK16Poly5.eps;
