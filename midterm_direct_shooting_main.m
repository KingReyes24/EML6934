%% EML6934 Optimal Control
%  Name:       Elias Reyes
%  Date:       04 April 2022
%  Assignment: Midterm
%  Goal:       Solve the orbit transfer problem using direct shooting
close all; clear all; clc
format longg
format compact

% specify number of coefficients to loop through
n_list = [ 2 3 4 5 6 ];

% specify number of figures per case
num_figs = 2;

% intitialize performance vectors
numn          = numel(n_list);
terminal_time = zeros(numn,1);
terminal_mass = zeros(numn,1);
iter          = zeros(numn,1);
sim_time      = zeros(numn,1);

% set number of states in problem
param.numStates = 5; % r,theta,vr,vtheta,m
% given conditions
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

count = 1; % counter for figures
% loop through different polynomial degrees cases
for idxn = 1:numn
    
    param.n         = n_list(idxn); % degree of polynomial used to parameterize the control
    param.numCoeff  = param.n+1;    % always 1 more coefficent then degree
    
    % initial guesses for unkown parameters 
    tfguess  = 10;
    cguess   = zeros(param.numCoeff,1);
    zguess   = [cguess; tfguess];

    % set up bounds
    cmin  = -50*ones(param.numCoeff,1);
    cmax  =  50*ones(param.numCoeff,1);
    tfmin = 0;
    tfmax = 100;
    zmin  = [cmin; tfmin];
    zmax  = [cmax; tfmax];
    A   = [];
    B   = [];
    Aeq = [];
    Beq = [];

    options = optimset('MaxFunEvals',100000,'MaxIter',1000,'Display','off');
    
    tic % start timer
    [solution,~,~,output] = fmincon(@orbitTransferObj,zguess,A,B,Aeq,Beq,zmin,zmax,@directOrbitTransferError,options,param);
    
    elapsed_time = toc;               % end timer and log
    iterations   = output.iterations; % log number of iterations

    %%  Take optimal initial conditions and integrate them for plotting purposes
    color   = ['b','r','g','c','m'];      
    tf      = solution(end);     % optimal tf
    c       = solution(1:end-1); % optimalcoefficients

    % initial conditions
    X0      = [param.r0; param.theta0; param.vr0; param.vtheta0; param.m0];

    % use tau scale [-1 1]
    tau     = [-1 1];
    options = odeset('reltol',1e-6);

    % propagate with the optimal conditions and coefficients to get the optimal trajectory 
    [t,p] = ode113(@directOrbitTransferOde,tau,X0,options,c,param,tf);

    % otpimal control
    beta  = polyval(c,t)*180/pi;    

    % plot the optimal trajectories and control
    figure(count); hold on; grid minor
    for p_idx = 1:param.numStates
        plot(t,p(:,p_idx),color(p_idx));
        h =  plot(t([1 end]),p([1 end],p_idx),'*k');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    xlabel('$\tau$','Interpreter','LaTeX')
    legend('$r(\tau)$','$\theta(\tau)$','$v_r(\tau)$','$v_\theta(\tau)$','$m(\tau)$','Interpreter','LaTeX')
    set(gcf,'color','white')
    set(gca,'fontweight','bold','fontsize',10)
    str1 = sprintf('States for trajectory that minimized fuel consumption');
    title(str1);

    figure(count +1); hold on; grid minor
    plot(t([1 end]),beta([1 end]),'*k')
    plot(t,beta)
    set(gcf,'color','white')
    set(gca,'fontweight','bold','fontsize',10)
    xlabel('$\tau$','Interpreter','LaTeX')
    ylabel('$\beta(\tau)$','Interpreter','LaTeX')
    str1 = sprintf('Optimal Control that minimizes fuel consumption');
    title(str1);
    
    % save off the data
    mf                  =  p(end,5);
    terminal_time(idxn) = tf;
    terminal_mass(idxn) = mf;
    iter(idxn)          = iterations;
    sim_time(idxn)      = elapsed_time;
    
    % incriment figure index
    count = count + num_figs;
end
