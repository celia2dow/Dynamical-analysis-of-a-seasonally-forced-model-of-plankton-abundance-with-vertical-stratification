% Script for producing the model figures in StAMPS article

% PLOT SETTINGS
clear
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
set(groot,'defaultFigureVisible','off')         % 'on' to turn back on
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
fontSize = 28;                                 % Font sizes with which to save each figure [21, 25, 28]
lwidth = 3;                                     % Width of plot lines
frame_width = 1.5;                              % Width of plot frame and ticks
dims_long = [0 0 42 22];                        % For setting all axes to have the same dimensions
dims_short = [0 0 25 22];
dims_smaller = [0 0 25 31];
dims_taller = [0 0 30 42];
legend_switch = 0;                              % 1: plot legend, 0: omit legend
tol = 0.001;                                     % Tolerance for the solution being steady in the long-term
% tolN = 0.09; tolP = 0.025; tolZ = 0.004;
tolN = 0.15; tolP = 0.03; tolZ = 0.01;            % Tolerance for jumps in the bifurcation diagrams of N, P and Z
tolLim = 0.00005;                                 % Tolerance for limiting cycle
rng(1)                                          % Seed the random number generator with 1
bi_frac = 0.95;                                 % Fraction of time-points discarded as burn-in       

% MODEL CHOICES: on == 1, off == 0
params.vert_strat = 0;
params.changing_MLD = 0;
params.num_reps = 100;                           % Number of years
params.num_days = params.num_reps*365.25;       % Number of days in year
params.days = 0:params.num_days-1;              % Array of days for solving

% PARAMETERS 
params.a = 0.184; % 1/(m day)                   % Related to the mazimum growth rate of P (0.184) (0.08)
% Alternative calculation of a:
% params.assumed_MLD = 12.5; % m                % Assumed MLD in non-chanding model
% params.Pmax = 2; % 1/(day) 0.2*12.5/2.58      % Maximum phytoplankton growth rate under optimal light conditions (0.5-2)
% params.a = 2.58*params.Pmax/params.assumed_MLD; % 1/(m day)
params.b = 0.23; % 1/m                          % Light attenuation by water (0.23) (0.1)
% Alternative calculation of b:
% params.assumed_Zeu = 100; % m                 % Assumed euphotic zone depth
% params.b = log(100)/params.assumed_Zeu; % 1/m
params.c = 0.4; % m^2/(g C)                     % P self-shading coefficient
params.d = 2.05; % m^3/(g C day)                % Higher predation on Z (for 1996 paper) (1.5 for steady state, 2.05 for oscillations)
params.e = 0.05; %0.05; % g C/(m^3)             % Half-saturation constant for N uptake
params.k = 0.05; % 1/day                        % Cross-thermocline exchange rate
params.r = 0.15; % 1/day                        % P respiration rate
params.s = 0.04; % 1/day                        % P sinking loss rate
params.N_0 = 0.8; % g C/(m^3)                   % N concentration below mixed layer
params.alpha = 0.25;                            % Z growth efficiency
params.beta = 0.33;                             % Z excretion fraction
params.gamma = 0.5;                             % Regeneration of Z predation excretion
params.lambda = 0.6; % 1/day                    % Maximum Z grazing rate
params.mu = 0.035; % g C/(m^3)                  % Z grazing half-saturation coefficient 

% SCALED BY ASSUMED MLD
params.assumed_MLD = 12.5; % m                  % Assumed MLD in non-changing model
params.assumed_fphiT = 1;  % °C m              % Assumed TSI in non-changing model
params.a_hat = params.a .* params.assumed_MLD;  % 1/m     
params.k_hat = params.k .* params.assumed_MLD;
params.s_hat = params.s .* params.assumed_MLD;

% TSI SCALING
params.TSI_mean = -66.7949; % °C m              % 10-year climatological mean temperature stratification index (TSI)
params.TSI_std_dev = 41.7331; % °C m            % Standard deviation of the above
params.strat_func = @(t) (1 + 200.^( ...         % Vertical stratification function given TSI at time t
    2.4 .* cos (2.*pi.*(t+90).*params.num_reps ./ params.num_days) ...
    - 0.4 .* params.TSI_mean / params.TSI_std_dev ...
    - 1)).^(-1);
% Alternative formulation of vertical stratification function:
% params.TSI_func = @(t) ...                    % Approximate value of TSI at time t given sinusoidal temperature 
%     1.2*(-params.TSI_std_dev*cos(2*pi*(t+90)/params.num_days)...
%     +params.TSI_mean);
% params.strat_func = @(t) ...                  % Vertical stratification function given TSI at time t
%     (1+exp(-log(40000)*(params.TSI_func(t)-params.TSI_mean)/params.TSI_std_dev ...
%     + log(1/sqrt(40000))))^(-1);



% FOLDER PATH
date_time = num2str(fix(clock));
folder_name = 'bifurcation_plots';
folder_path = [pwd '/' date '/' folder_name];
% The folder paths used agrees with a macOS system and may need to be
% flipped to '\' for a Windows system
if ~exist(folder_path, 'dir')
    mkdir(folder_path)
end

% COLOURS: uisetcolor
col.grey = [0.6510 0.6510 0.6510];
col.nutri_orange = [1.0000 0.4118 0.1608];
col.phyto_green = [0.0157 0.7804 0.4235]; 
col.zoop_blue = [0.0275 0.4471 0.9294];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME-SERIES TRAJECTORIES
% PREPARE INITIAL CONDITION
N0s = 0.4; % N in [0.2, 0.4] g C / m^3
P0s = 0.1; % P in [0.05, 0.1] g C / m^3
Z0s = 0.05; % Z in [0.05, 0.055] g C / m^3

% SOLVE FOR ODE CONCENTRATIONS
[t,y] = ode45(@(t,x) odefunc(t,x,params),...
    [params.days(1) params.days(end)],y0); % or ode15s, optionally with options

times=t; 
Narray = y(:, 1); 
Parray = y(:, 2);
Zarray = y(:, 3);

if params.num_reps>1
    timesTT = times(times>=(params.num_reps-1)*365.25)-(params.num_reps-1)*365;
    NarrayTT = Narray(times>=(params.num_reps-1)*365.25);
    ParrayTT = Parray(times>=(params.num_reps-1)*365.25);
    ZarrayTT = Zarray(times>=(params.num_reps-1)*365.25);
else
    timesTT = times;
    NarrayTT = Narray;
    ParrayTT = Parray;
    ZarrayTT = Zarray;
end

close all

% PLOT TIME-TRAJECTORIES
fig1 = figure(1);

plot(timesTT, NarrayTT, 'Color', col.nutri_orange, 'LineWidth', lwidth); % Nutrient
hold on
plot(timesTT, ParrayTT, 'Color', col.phyto_green, 'LineWidth', lwidth);  % Phytoplankton
plot(timesTT, ZarrayTT, 'Color', col.zoop_blue, 'LineWidth', lwidth);    % Zooplankton

ylabel("Concentration (g C/m$^3$)")
xlabel("Time (days)")
xlim([0,timesTT(end)])
if legend_switch
    leg1 = legend('Mixed layer nutrient ($N$)', 'Phytoplankton ($P$)', 'Zooplankton ($Z$)');
    leg1.Location = 'east';
    leg1.Box ="off";
end
hold off

fontsize(fontSize,'points')
set(fig1, 'Units', 'centimeters', 'OuterPosition', dims_long)
exportgraphics(fig1, folder_path + "/time_trajectory_d" + ...
    num2str(params.d) + "_r" + num2str(params.r) + ...
    "_vertStrat" + num2str(params.vert_strat) + ...
    "_changeMLD" + num2str(params.changing_MLD) + ...
    ".pdf",'ContentType','vector')

% PLOT PHASE PORTRAIT
fig2 = figure(2);

% Find burn in
total_tsteps = length(times);
if params.num_reps>1
    burn_in = find(times>=(params.num_reps-1)*365.35,1);
elseif abs(Parray(end)-Parray(end-1))<tol
    burn_in = floor(3*total_tsteps/4);  
else
    burn_in = floor(total_tsteps/9);
end

% Plot burn in
plot3(Parray(1:burn_in), Narray(1:burn_in), Zarray(1:burn_in), ...
    'Color', [col.grey, 0.5], 'LineWidth', lwidth)
hold on

% Plot long-time behaviour
plot3(Parray(burn_in:end), Narray(burn_in:end), Zarray(burn_in:end), ...
    'Color', 'k', 'LineWidth', lwidth)
hold off

xlim([0,max(Parray,[],"all")])
ylim([0,max(Narray,[],"all")])
zlim([0,max(0.2, max(Zarray,[],"all"))])
grid on
ylabel('$N$ (g C/m$^3$)')
xlabel('$P$ (g C/m$^3$)')
zlabel('$Z$ (g C/m$^3$)')
set(gca, 'YDir','reverse')
if legend_switch
    leg2 = legend('Transitory behaviour', 'Long-term behaviour');
    leg2.Location = 'north';
    leg2.Box ="off";
end

fontsize(fontSize,'points')
set(fig2, 'Units', 'centimeters', 'OuterPosition', dims_short)
exportgraphics(fig2,folder_path + "/phase_plot_d" + ...
    num2str(params.d) + "_r" + num2str(params.r) + ...
    "_vertStrat" + num2str(params.vert_strat) + ...
    "_changeMLD" + num2str(params.changing_MLD) + ...
    ".pdf",'ContentType','vector')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIFURCATION DIAGRAMS IN r 
% PREPARE INITIAL CONDITION
num_ICs = 30;
N0s = makeSeqVec(0.001,1,num_ICs);     % [0.05,0.45] N g C / m^3
P0s = makeSeqVec(0.001,1,num_ICs);   % [0.025,0.325] P g C / m^3
Z0s = makeSeqVec(0.001,0.2,num_ICs);     % [0.03,0.08] Z g C / m^3

% MAKE BIRUFCATION PLOT r and Phi
params.r_vec = 0:0.005:0.8; %0.06:0.001:0.15;
params.sHat_vec = 0:0.005:(0.08*12.5);
params.kHat_vec = 0:0.005:(0.13*12.5);
params.aHat_vec = 0:0.005:(0.28*12.5);

% ITERATE THROUGH BIFURCATION PARAMETER VALUES AND FIND LONG TERM BEHAVIOUR
param_i=0; % Initialise index for array of bifurcation paremeter values
vec = params.r_vec; % CHANGE DEPENDING ON WHICH BIFURCATION PLOT
len_vec = length(vec);
params.Phi_vec =(params.a*params.N_0/(params.b*(params.e+params.N_0))-params.s-params.k)-vec;

% Initialise storage for all unordered solutions
Nstore = nan(len_vec,num_ICs,6); % Col 1: fsolve solution, Col 2: fsolve solution stability,
Pstore = nan(len_vec,num_ICs,6); % Col 3: ode45 solution minimum, 4: ode45 solution maximum,
Zstore = nan(len_vec,num_ICs,6); % Col 5: ode45 limit existence, Col 6: ode45 solution stability,

% Initialise storage for sorted stable and unstable points and limit cycles
% Each field will contain a list of points and their stabilities
Nstore_wStability.stable.sol_1 = []; Nstore_wStability.unstable.sol_1 = []; 
Nstore_wStability.stable.sol_ss2 = []; Nstore_wStability.unstable.sol_ss2 = []; 
Nstore_wStability.stable.sol_lc2 = []; Nstore_wStability.unstable.sol_lc2 = []; 

Pstore_wStability.stable.sol_1 = []; Pstore_wStability.unstable.sol_1 = []; 
Pstore_wStability.stable.sol_ss2 = []; Pstore_wStability.unstable.sol_ss2 = []; 
Pstore_wStability.stable.sol_lc2 = []; Pstore_wStability.unstable.sol_lc2 = []; 

Zstore_wStability.stable.sol_1 = []; Zstore_wStability.unstable.sol_1 = []; 
Zstore_wStability.stable.sol_ss2 = []; Zstore_wStability.unstable.sol_ss2 = []; 
Zstore_wStability.stable.sol_lc2 = []; Zstore_wStability.unstable.sol_lc2 = []; 

% ITERATE THROUGH BIFURCATION PARAMETER VALUES
for param_i = 1:len_vec
    % For tracking
    fprintf('%d\n',param_i)

    % Set current value of bifurcation parameter
    params.r = vec(param_i);  % CHANGE DEPENDING ON WHICH BIFURCATION PLOT

    % ITERATE THROUGH INITIAL CONDITIONS
    N0s_parami = N0s; P0s_parami = P0s; Z0s_parami = Z0s;
    for ic = 1:num_ICs
        % Pick N0 and remove it from the list
        N0pick = randi([1 num_ICs - ic + 1],1); 
        N0 = N0s_parami(N0pick); 
        N0s_parami(N0pick)=[]; 

        % Pick P0 and remove it from the list
        P0pick = randi([1 num_ICs - ic + 1],1); 
        P0 = P0s_parami(P0pick); 
        P0s_parami(P0pick)=[]; 

        % Pick Z0 and remove it from the list
        Z0pick = randi([1 num_ICs - ic + 1],1); 
        Z0 = Z0s_parami(Z0pick); 
        Z0s_parami(Z0pick)=[]; 

        y0 = [N0, P0, Z0];

        % FSOLVE: Solve the system of ODEs for a steady state from given IC
        odefunc_time1 = @(x) odefunc(1,x,params);
        options = optimoptions('fsolve','MaxIterations',300,...
            'StepTolerance',1e-12,'Display','none'); % Reduce the number of permitted iterations
        [x_longterm,fval,exitflag,output,jacobian]=...
            fsolve(odefunc_time1,y0,options);

        % Abort iteration if ...
        if exitflag ~=1 ...                         % ... fsolve did not converge,
                || any(x_longterm <= 0) ...         % ... the solution is negative,
                || any(x_longterm >= params.N_0)    % ... or the solution is too large.
            continue
        end

        % Otherwise, continue through current iteration
        Nstore(param_i,ic,1) = x_longterm(1);
        Pstore(param_i,ic,1) = x_longterm(2);
        Zstore(param_i,ic,1) = x_longterm(3);

        % Check if FSOLVE solution is stable on not
        if any(real(eig(jacobian))>=0)  % Solution is unstable
            Nstore(param_i,ic,2) = 0;    
            Pstore(param_i,ic,2) = 0;
            Zstore(param_i,ic,2) = 0;
        else                            % Solution is stable
            Nstore(param_i,ic,2) = 1;   % Indicate steady state solution is stable
            Pstore(param_i,ic,2) = 1;
            Zstore(param_i,ic,2) = 1;
        end  

        % ODE45: Solve the system of ODEs for time series solutions from given IC
        [~,y] = ode45(@(t,x) odefunc(t,x,params),...
                [params.days(1) params.days(end)],y0); % or ode15s, optionally with options
        [minN, maxN, isLimitN, Nstore] = findFinalPeaks(y, 1, bi_frac, tolLim, Nstore, param_i, ic);
        [minP, maxP, isLimitP, Pstore] = findFinalPeaks(y, 2, bi_frac, tolLim, Pstore, param_i, ic);
        [minZ, maxZ, isLimitZ, Zstore] = findFinalPeaks(y, 3, bi_frac, tolLim, Zstore, param_i, ic);
    end

    % DETERMINE SOLUTION STABILITY
    Nstore_wStability = collectSolsPerBifurcVal(Nstore, param_i, tol, vec, Nstore_wStability,1);
    Pstore_wStability = collectSolsPerBifurcVal(Pstore, param_i, tol, vec, Pstore_wStability,2);
    Zstore_wStability = collectSolsPerBifurcVal(Zstore, param_i, tol, vec, Zstore_wStability,3);
end

% PLOT BIFURCATION DIAGRAMS
close all
Nruns = makeAndSaveBifurcationPlot(Nstore_wStability, 1, tolN, vec, ...
    col.nutri_orange, lwidth, fontSize, params, dims_smaller, folder_path);
Pruns = makeAndSaveBifurcationPlot(Pstore_wStability, 2, tolP, vec, ...
    col.phyto_green, lwidth, fontSize, params, dims_smaller, folder_path);
Zruns = makeAndSaveBifurcationPlot(Zstore_wStability, 3, tolZ, vec, ...
    col.zoop_blue, lwidth, fontSize, params, dims_smaller, folder_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIFURCATION DIAGRRAM SEASONALLY FORCED
% PREPARE INITIAL CONDITION
num_ICs = 30;
N0s = makeSeqVec(0.001,1,num_ICs);     % [0.05,0.45] N g C / m^3
P0s = makeSeqVec(0.001,1,num_ICs);   % [0.025,0.325] P g C / m^3
Z0s = makeSeqVec(0.001,0.2,num_ICs);     % [0.03,0.08] Z g C / m^3

% MAKE BIRUFCATION PLOT - SEASONALLY FORCED
params.times = 0:0.75:365.25;
len_vec = length(params.times);
params.strat_vec = params.strat_func(params.times);
params.M_vec = zeros(1,len_vec);
params.h_vec = zeros(1,len_vec);
for i=1:length(params.times)
    [depth, h] = MLD_func(params.times(i));
    params.M_vec(i) = depth;    
    params.h_vec(i) = h;
end
params.A_vec = params.a_hat ./ params.M_vec;
params.S_vec = params.s_hat .* params.strat_vec ./ params.M_vec;
params.K_vec = params.k_hat .* params.strat_vec ./ params.M_vec;

% ITERATE THROUGH BIFURCATION PARAMETER VALUES AND FIND LONG TERM BEHAVIOUR
params.PhiNew_vec = params.A_vec*params.N_0/(params.b*(params.e+params.N_0))-params.r-params.S_vec-params.K_vec - max(0,params.h_vec)./params.M_vec;
vec = params.times;

% Initialise storage for all unordered solutions
Nstore = nan(len_vec,num_ICs,6); % Col 1: fsolve solution, Col 2: fsolve solution stability,
Pstore = nan(len_vec,num_ICs,6); % Col 3: ode45 solution minimum, 4: ode45 solution maximum,
Zstore = nan(len_vec,num_ICs,6); % Col 5: ode45 limit existence, Col 6: ode45 solution stability,

% Initialise storage for sorted stable and unstable points and limit cycles
% Each field will contain a list of points and their stabilities
Nstore_wStability.stable.sol_1 = []; Nstore_wStability.unstable.sol_1 = []; 
Nstore_wStability.stable.sol_ss2 = []; Nstore_wStability.unstable.sol_ss2 = []; 
Nstore_wStability.stable.sol_lc2 = []; Nstore_wStability.unstable.sol_lc2 = []; 

Pstore_wStability.stable.sol_1 = []; Pstore_wStability.unstable.sol_1 = []; 
Pstore_wStability.stable.sol_ss2 = []; Pstore_wStability.unstable.sol_ss2 = []; 
Pstore_wStability.stable.sol_lc2 = []; Pstore_wStability.unstable.sol_lc2 = []; 

Zstore_wStability.stable.sol_1 = []; Zstore_wStability.unstable.sol_1 = []; 
Zstore_wStability.stable.sol_ss2 = []; Zstore_wStability.unstable.sol_ss2 = []; 
Zstore_wStability.stable.sol_lc2 = []; Zstore_wStability.unstable.sol_lc2 = []; 

% ITERATE THROUGH BIFURCATION PARAMETER VALUES
for param_i = 1:len_vec
    % For tracking
    fprintf('%d\n',param_i)

    % Set current value of seasonally forced parameters
    params.assumed_MLD = params.M_vec(param_i);
    params.assumed_fphiT = params.strat_vec(param_i);

    % ITERATE THROUGH INITIAL CONDITIONS
    N0s_parami = N0s; P0s_parami = P0s; Z0s_parami = Z0s;
    for ic = 1:num_ICs
        % Pick N0 and remove it from the list
        N0pick = randi([1 num_ICs - ic + 1],1); 
        N0 = N0s_parami(N0pick); 
        N0s_parami(N0pick)=[]; 

        % Pick P0 and remove it from the list
        P0pick = randi([1 num_ICs - ic + 1],1); 
        P0 = P0s_parami(P0pick); 
        P0s_parami(P0pick)=[]; 

        % Pick Z0 and remove it from the list
        Z0pick = randi([1 num_ICs - ic + 1],1); 
        Z0 = Z0s_parami(Z0pick); 
        Z0s_parami(Z0pick)=[]; 
        
        y0 = [N0, P0, Z0];
        
        % FSOLVE: Solve the system of ODEs for a steady state from given IC
        odefunc_time1 = @(x) odefunc(1,x,params);
        options = optimoptions('fsolve','MaxIterations',300,...
            'StepTolerance',1e-12,'Display','none'); % Reduce the number of permitted iterations
        [x_longterm,fval,exitflag,output,jacobian]=...
            fsolve(odefunc_time1,y0,options);
        
        % Abort iteration if ...
        if exitflag ~=1 ...                         % ... fsolve did not converge,
                || any(x_longterm <= 0) ...         % ... the solution is negative,
                || any(x_longterm >= params.N_0)    % ... or the solution is too large.
            continue
        end

        % Otherwise, continue through current iteration
        Nstore(param_i,ic,1) = x_longterm(1);
        Pstore(param_i,ic,1) = x_longterm(2);
        Zstore(param_i,ic,1) = x_longterm(3);

        % Check if FSOLVE solution is stable on not
        if any(real(eig(jacobian))>=0)  % Solution is unstable
            Nstore(param_i,ic,2) = 0;    
            Pstore(param_i,ic,2) = 0;
            Zstore(param_i,ic,2) = 0;
        else                            % Solution is stable
            Nstore(param_i,ic,2) = 1;   % Indicate steady state solution is stable
            Pstore(param_i,ic,2) = 1;
            Zstore(param_i,ic,2) = 1;
        end  

        % ODE45: Solve the system of ODEs for time series solutions from given IC
        [~,y] = ode45(@(t,x) odefunc(t,x,params),...
                [params.days(1) params.days(end)],y0); % or ode15s, optionally with options
        [minN, maxN, isLimitN, Nstore] = findFinalPeaks(y, 1, bi_frac, tolLim, Nstore, param_i, ic);
        [minP, maxP, isLimitP, Pstore] = findFinalPeaks(y, 2, bi_frac, tolLim, Pstore, param_i, ic);
        [minZ, maxZ, isLimitZ, Zstore] = findFinalPeaks(y, 3, bi_frac, tolLim, Zstore, param_i, ic);
    end

    % DETERMINE SOLUTION STABILITY
    Nstore_wStability = collectSolsPerBifurcVal(Nstore, param_i, tol, vec, Nstore_wStability);
    Pstore_wStability = collectSolsPerBifurcVal(Pstore, param_i, tol, vec, Pstore_wStability);
    Zstore_wStability = collectSolsPerBifurcVal(Zstore, param_i, tol, vec, Zstore_wStability);
end

% PLOT BIFURCATION DIAGRAMS
close all

Nruns = makeAndSaveBifurcationPlot(Nstore_wStability, 1, tolN, vec, ...
    col.nutri_orange, lwidth, fontSize, params, dims_smaller, folder_path);
Pruns = makeAndSaveBifurcationPlot(Pstore_wStability, 2, tolP, vec, ...
    col.phyto_green, lwidth, fontSize, params, dims_smaller, folder_path);
Zruns = makeAndSaveBifurcationPlot(Zstore_wStability, 3, tolZ, vec, ...
    col.zoop_blue, lwidth, fontSize, params, dims_smaller, folder_path);