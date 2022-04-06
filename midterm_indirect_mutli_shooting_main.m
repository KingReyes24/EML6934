%% EML6934 Optimal Control
%  Name:       Elias Reyes
%  Date:       04 April 2022
%  Assignment: Midterm
%  Goal:       Solve the orbit transfer problem using indirect multiple shooting
close all; clear all; clc
format longg
format compact

% specify number of intervals to loop through
k_list = [ 2 4 8 16];

% specify number of figures per case
num_figs = 2;

% intitialize performance vectors
numk          = numel(k_list);
terminal_time = zeros(numk,1);
terminal_mass = zeros(numk,1);
iter          = zeros(numk,1);
sim_time      = zeros(numk,1);

% set num states
param.numStates = 10;

% given constants
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

% transversality conditions
param.lamthetaf = 0;
param.lammf     = 1;
param.Hf        = 0; 

% initial guesses for unkown co-states 
lamr            = -2;
lamtheta        = 0;
lamvtheta       = 2;
lamvr           = 2;
lamm            = 2;
tf              = 3.5;

count         = 1; % count for figs
% loop through all the cases
for kidx = 1:numel(k_list)
    % set specific interval case
    param.k         = k_list(kidx);
    
    % Create conditions that need to be solved. This takes the intial guesses
    % above and repeats them for the desired number of intervals. This also sets
    % the guess for the states at each interval. 
    lamguess   = [lamr; lamtheta; lamvtheta; lamvr; lamm];
    stateguess = ones(param.numStates/2,1);
    pguess     = repmat([stateguess;lamguess],param.k-1,1); % state/cp-state guesses
    tfguess    = tf;
    % tf guess
    % concatinate all the unkowns. Begin with the only the costates since we know
    % the states at t0
    zguess     = [lamguess; pguess; tfguess];               

    options = optimoptions('fsolve','Display','off','TolFun',1e-8);
    
    % solve the optimal control problem
    f = @(x)indirectMultiOrbitTransferError(x,param);

    tic % start timer
    [solution,~,~,output] = fsolve(f, zguess ,options);

    elapsed_time = toc; % end timer and log
    iterations   = output.iterations;

    %%  Take optimal initial conditions and integrate them for plotting purposes
    color   = ['b','r','g','c','m'];
    r0      = param.r0;
    theta0  = param.theta0;
    vr0     = param.vr0;
    vtheta0 = param.vtheta0;
    m0      = param.m0;

    tf      = solution(end); % optimal tf
    % optimal conditions at each interval
    zvec    = [r0; theta0; vr0; vtheta0; m0; solution(1:end-1)];
    % reshape to make number of columns equal to number of intervals. This
    % helps with indexing.
    zvec    = reshape(zvec,param.numStates,param.k);
    time    = linspace(param.t0,tf,param.k+1);
    options = odeset('reltol',1e-6);

    for idx = 1:param.k
        X0    = zvec(:,idx);
        tspan = [time(idx) time(idx+1)];

        [t,p] = ode113(@orbitTransferOde,tspan,X0,options,param);

        s = length(t);

        lamvrf     = p(:,8);
        lamvthetaf = p(:,9);

        % optimal control 
        beta = mod(atan2(lamvrf,lamvthetaf),2*pi);

        % plot figures
        figure(count); hold on; 
        for p_idx = 1:param.numStates/2
            plot(t,p(:,p_idx),color(p_idx));
            h =  plot(t([1 end]),p([1 end],p_idx),'*k');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        xlabel('t','Interpreter','LaTeX')
        legend('r(t)','$\theta(t)$','vr(t)','v$\theta(t)$','m(t)','Interpreter','LaTeX')
        set(gcf,'color','white')
        set(gca,'fontweight','bold','fontsize',10,'XMinorGrid','on','YMinorGrid','on')
        str1 = sprintf('States for trajectory that minimized fuel consumption (%d Intervals)',param.k);
        title(str1);

        figure(count+1); hold on;
        plot(t([1 end]),beta([1 end])*180/pi,'*k')
        plot(t,beta*180/pi)
        set(gcf,'color','white')
        set(gca,'fontweight','bold','fontsize',10,'XMinorGrid','on','YMinorGrid','on')
        xlabel('t','Interpreter','LaTeX')
        ylabel('$\beta(t)$','Interpreter','LaTeX')
        str1 = sprintf('Optimal Control that minimizes fuel consumption (%d Intervals)',param.k);
        title(str1);  
    end
    % save off the data
    mf                  =  p(end,5);
    terminal_time(kidx) = tf;
    terminal_mass(kidx) = mf;
    iter(kidx)          = iterations;
    sim_time(kidx)      = elapsed_time;
    
    % incriment figure index
    count = count + num_figs;
end

