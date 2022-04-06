%% EML6934 Optimal Control
%  Name:       Elias Reyes
%  Date:       04 April 2022
%  Assignment: Midterm
%  Goal:       Solve the orbit transfer problem using direct multiple shooting
close all; clear all; clc
format longg
format compact

% specify desired polynomial degrees and intervals
n_list = [2 3 ] %4 5 6];
k_list = [2 4 ] %8 16];

% specify number of figures in loop
num_figs = 2;

% initialize performance matrices
numk          = numel(k_list);
numn          = numel(n_list);
terminal_time = zeros(numn,numk);
terminal_mass = zeros(numn,numk);
iter          = zeros(numn,numk);
sim_time      = zeros(numn,numk);

% set number of states
param.numStates = 5; % r,theta,vr,vtheta,m 

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

% intitialize counter for plots
count           = 1;

for idxn = 1:numn
    
    % degree of polynomial used to parameterize the control
    param.n         = n_list(idxn); 
    param.numCoeff  = param.n+1; % always 1 more coefficent then degree
    
    for idxk = 1:numk
        
        % number of intervals
        param.k         = k_list(idxk); 

        % initial guesses for unkown parameters 
        tfguess  = 10;
        pguess   = ones(param.numStates*(param.k-1),1);
        cguess   = zeros(param.numCoeff*param.k,1);
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

        options = optimset('MaxFunEvals',100000,'MaxIter',1000,'Display','off','TolFun',1e-4);

        tic % starter timer
        [solution,~,~,output] = fmincon(@orbitTransferObj,zguess,A,B,Aeq,Beq,zmin,zmax,@directMultiOrbitTransferError,options,param);
        
        elapsed_time = toc; % stop timer and record
        iterations   = output.iterations; % record num iterations

        %%  Take optimal conditions and integrate them for plotting purposes
        color   = ['b','r','g','c','m'];
        % states at t0
        r0      = param.r0;
        theta0  = param.theta0;
        vr0     = param.vr0;
        vtheta0 = param.vtheta0;
        m0      = param.m0;
        P0      = [r0; theta0; vr0; vtheta0; m0];

        % seperate the unknowns 
        % P is the unknown states
        P_end  = param.numStates*(param.k-1);
        P_tmp  = solution(1:P_end);
        P      = [P0;P_tmp];
        P      = reshape(P,param.numStates,[]); % reshape to make each column an interval
        % c is the uknown coeffients
        c_end  = numel(solution)-1;
        c_list = solution(P_end+1:c_end);
        c_list = reshape(c_list,param.numCoeff,[]); % reshape to make each column an interval
        tf     = solution(end); 

        % create tau grid
        tau     = linspace(-1,1,param.k+1);
        options = odeset('reltol',1e-6);

        for idx = 1:param.k
            % extract coefficents for the specific interval
            c  = c_list(:,idx); 
            X0 = P(:,idx);
            % tau span for specific interval
            tspan = [tau(idx) tau(idx+1)];
            [t,p] = ode113(@directOrbitTransferOde,tspan,X0,options,c,param,tf);

            % optimal control
            beta =  mod(polyval(c,t),2*pi);

            % plot figures
            figure(count); hold on;
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
            
            figure(count+1); hold on; 
            plot(t([1 end]),beta([1 end])*180/pi,'*k')
            plot(t,beta*180/pi)
            set(gcf,'color','white')
            set(gca,'fontweight','bold','fontsize',10,'XMinorGrid','on','YMinorGrid','on')
            xlabel('$\tau$','Interpreter','LaTeX')
            ylabel('$\beta(\tau)$','Interpreter','LaTeX')
            str1 = sprintf('Optimal Control that minimizes fuel consumption (%d Intervals)',param.k);
            title(str1);
        end
        
        % save off data
        mf                       =  p(end,5);
        terminal_time(idxn,idxk) = tf;
        terminal_mass(idxn,idxk) = mf;
        iter(idxn,idxk)          = iterations;
        sim_time(idxn,idxk)      = elapsed_time;
        
        % incriment figure index
        count = count + num_figs;
    end
end

%% Put all the data into a stuct 
fields = {'N','K','Iterations','SimTime','TerminalTime','TerminalMass'};
count  = 1;
for idxn = 1:numn
    N = n_list(idxn);
    for idxk = 1:numk
        fieldname = sprintf('combo%d',count);
        K   = k_list(idxk);
        itr = iter(idxn,idxk);
        t   = sim_time(idxn,idxk);
        tf  = terminal_time(idxn,idxk);
        mf  = terminal_mass(idxn,idxk);
        stuff = [N K itr t tf mf];
        for kfield = 1:numel(fields)
            junk = stuff(kfield);
            data.(fieldname).(fields{kfield}) = junk;
        end
        count = count + 1;
    end
end

%% Additional Figures
nlistvec = repmat(n_list',1,numel(k_list));
klistvec = repmat(k_list,numel(n_list),1);

figure;
surf(nlistvec,klistvec,sim_time, 'FaceAlpha',0.75)
xlabel('Polynomial Degree N')
ylabel('Numer of Intervals K')
zlabel('Execution Time')
set(gcf,'color','white')
set(gca,'fontweight','bold','fontsize',10, 'XMinorGrid','on','YMinorGrid','on')
title('Direct Multiple Shooting Performance')

figure; hold on;
for idx = 1:numn
plot(klistvec(idx,:),sim_time(idx,:),'-o')
end
set(gcf,'color','white')
set(gca,'fontweight','bold','fontsize',10, 'XMinorGrid','on','YMinorGrid','on')
xlabel('Intervals')
ylabel('Execution Time')
legend('Degree 2','Degree 3','Degree 4','Degree 5','Degree 6')
title('Execution Performance for Direct Multi Shooting')