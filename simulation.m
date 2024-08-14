% Run model simulations for Archibald et al. (2024), The American Naturalist
% 9.9.0.1538559 (R2020b) Update 3

%% Simulation 1: Reciprocal transplant simulations (manuscript Fig. 2)

% Assign save directory
sd1 = 'fig2_data';

% Set up simulation
tspan = 1:0.1:20; % time span
I = 3.0*ones(length(tspan),1); % light
T = 16:1:32; % temperature

% Load model parameters file
p = params;
p(13) = 0; % turn off evolution
Tevo = 16:1:32; % set range of evolved temperatures

% Set up optimal strategy update framework
tau = 10; % update frequency
opt_time = min(tspan):(1/tau):max(tspan);
global opt_R

% ODE solver options
opt = odeset('RelTol', 1e-4, 'AbsTol', 1e-8, 'NonNegative', 1:8);

% Run simulations and calculate diagnostics
dn = cell(length(T),length(Tevo));
for j = 1:length(Tevo)
    x0 = [1e6 0.15 0.15 0.7 Tevo(j) 10 1.25 1e3]; % initial conditions
    for i = 1:length(T)
        opt_R = NaN(length(opt_time),3);
        [t, x] = ode45(@(t,y) mixo_model(t, y, tspan, T(i)*ones(length(tspan),1), I, p, opt_time), tspan, x0, opt);
        dn{i,j} = mod_diag(x,t,p,T(i),I(1));
    end
end
save(sd1, 'T', 'Tevo', 'dn');

%% Simulation 2: Tempertaure increase control timeseries (manuscript Fig. 3)

% Assign save directory
sd2 = 'fig3_data_tinc_con';

% Set up simulation
tspan = 0:0.1:1; % dilution period
n = 200; % length of simulation (days)
I = 3.0*ones(n*length(tspan),1); % light
T = [24*ones(floor(0.25*length(I)),1); 30*ones(ceil(0.75*length(I)),1)]; % temperature

% Load parameters
p = params;
p(13) = 0; % turn off evolution

% Set up optimal strategy update framework
tau = 10;
opt_time = min(tspan):(1/tau):max(tspan);
global opt_R

% ODE solver options
opt = odeset('RelTol', 1e-4, 'AbsTol', 1e-8, 'NonNegative', 1:7);

% Run simulation with evolution
x0 = [1e6 0.15 0.15 0.7 24 10 1.25 1e3]; % initial conditions
x = NaN(n*length(tspan),length(x0));
for i = 1:n
    if i > 1
        x0 = xtemp(end,:);
        % Dilution criteria
        if x0(8) > 1e6
            bm = x0(1)/x0(8);
            x0(1) = 1e8;
            x0(8) = 1e3;
        end
    end
    ii = (i-1)*length(tspan)+1:i*length(tspan);
    opt_R = NaN(length(opt_time),3);
    [ttemp, xtemp] = ode45(@(t,y) mixo_model(t, y, tspan, T(ii), I(ii), p, opt_time), tspan, x0, opt);
    if i == 1
        t = ttemp;
    else
        t(ii) = t(end) + ttemp;
    end
    x(ii,:) = xtemp;
end
dn_tinc_con = ts_diag(x,t,p,T,I);
Tinc = T;
save(sd2, 'Tinc', 'dn_tinc_con');

%% Simulation 3: Tempertaure increase evolving timeseries (manuscript Fig. 3)

% Assign save directory
sd3 = 'fig3_data_tinc_evo';

% Set up simulation
tspan = 0:0.1:1; % dilution period
n = 200; % length of simulation (days)
I = 3.0*ones(n*length(tspan),1); % light
T = [24*ones(floor(0.25*length(I)),1); 30*ones(ceil(0.75*length(I)),1)]; % temperature

% Load parameters
p = params;

% Set up optimal strategy update framework
tau = 10;
opt_time = min(tspan):(1/tau):max(tspan);
global opt_R

% ODE solver options
opt = odeset('RelTol', 1e-4, 'AbsTol', 1e-8, 'NonNegative', 1:7);

% Run simulation with evolution
x0 = [1e6 0.15 0.15 0.7 24 10 1.25 1e3]; % initial conditions
x = NaN(n*length(tspan),length(x0));
for i = 1:n
    if i > 1
        x0 = xtemp(end,:);
        % Dilution criteria
        if x0(8) > 1e6
            bm = x0(1)/x0(8);
            x0(1) = 1e8;
            x0(8) = 1e3;
        end
    end
    ii = (i-1)*length(tspan)+1:i*length(tspan);
    opt_R = NaN(length(opt_time),3);
    [ttemp, xtemp] = ode45(@(t,y) mixo_model(t, y, tspan, T(ii), I(ii), p, opt_time), tspan, x0, opt);
    if i == 1
        t = ttemp;
    else
        t(ii) = t(end) + ttemp;
    end
    x(ii,:) = xtemp;
end
dn_tinc_evo = ts_diag(x,t,p,T,I);
save(sd3, 'dn_tinc_evo');

%% Simulation 4: Tempertaure decrease control timeseries (manuscript Fig. 3)

% Assign save directory
sd4 = 'fig3_data_tdec_con';

% Set up simulation
tspan = 0:0.1:1; % dilution period
n = 200; % length of simulation (days)
I = 3.0*ones(n*length(tspan),1); % light
T = [24*ones(floor(0.25*length(I)),1); 18*ones(ceil(0.75*length(I)),1)]; % temperature

% Load parameters
p = params;
p(13) = 0; % turn off evolution

% Set up optimal strategy update framework
tau = 10;
opt_time = min(tspan):(1/tau):max(tspan);
global opt_R

% ODE solver options
opt = odeset('RelTol', 1e-4, 'AbsTol', 1e-8, 'NonNegative', 1:7);

% Run simulation with evolution
x0 = [1e6 0.15 0.15 0.7 24 10 1.25 1e3]; % initial conditions
x = NaN(n*length(tspan),length(x0));
for i = 1:n
    if i > 1
        x0 = xtemp(end,:);
        % Dilution criteria
        if x0(8) > 1e6
            bm = x0(1)/x0(8);
            x0(1) = 1e8;
            x0(8) = 1e3;
        end
    end
    ii = (i-1)*length(tspan)+1:i*length(tspan);
    opt_R = NaN(length(opt_time),3);
    [ttemp, xtemp] = ode45(@(t,y) mixo_model(t, y, tspan, T(ii), I(ii), p, opt_time), tspan, x0, opt);
    if i == 1
        t = ttemp;
    else
        t(ii) = t(end) + ttemp;
    end
    x(ii,:) = xtemp;
end
dn_tdec_con = ts_diag(x,t,p,T,I);
Tdec = T;
save(sd4, 'Tdec', 'dn_tdec_con');

%% Simulation 5: Tempertaure decrease evolving timeseries (manuscript Fig. 3)

% Assign save directory
sd5 = 'fig3_data_tdec_evo';

% Set up simulation
tspan = 0:0.1:1; % dilution period
n = 200; % length of simulation (days)
I = 3.0*ones(n*length(tspan),1); % light
T = [24*ones(floor(0.25*length(I)),1); 18*ones(ceil(0.75*length(I)),1)]; % temperature

% Load parameters
p = params;

% Set up optimal strategy update framework
tau = 10;
opt_time = min(tspan):(1/tau):max(tspan);
global opt_R

% ODE solver options
opt = odeset('RelTol', 1e-4, 'AbsTol', 1e-8, 'NonNegative', 1:7);

% Run simulation with evolution
x0 = [1e6 0.15 0.15 0.7 24 10 1.25 1e3]; % initial conditions
x = NaN(n*length(tspan),length(x0));
for i = 1:n
    if i > 1
        x0 = xtemp(end,:);
        % Dilution criteria
        if x0(8) > 1e6
            bm = x0(1)/x0(8);
            x0(1) = 1e8;
            x0(8) = 1e3;
        end
    end
    ii = (i-1)*length(tspan)+1:i*length(tspan);
    opt_R = NaN(length(opt_time),3);
    [ttemp, xtemp] = ode45(@(t,y) mixo_model(t, y, tspan, T(ii), I(ii), p, opt_time), tspan, x0, opt);
    if i == 1
        t = ttemp;
    else
        t(ii) = t(end) + ttemp;
    end
    x(ii,:) = xtemp;
end
dn_tdec_evo = ts_diag(x,t,p,T,I);
save(sd5, 'dn_tdec_evo');

%% Simulation 6: Sensitivity Analysis (supplement Fig. S1)

% Assign save directory
sd6 = 'figs1';

% Setup simulation
tspan = 1:0.1:10;
tau = 10;
I = 3.0*ones(length(tspan),1);
T = 18:2:30;
That = 18:2:30;
opt_time = min(tspan):(1/tau):max(tspan);
global opt_R
opt = odeset('RelTol', 1e-4, 'AbsTol', 1e-8, 'NonNegative', 1:7);

% Baseline parameters
u = NaN(length(T),length(That));
p = params;
p(13) = 0; % turn off evolution
for j = 1:length(That)
    x0 = [1e6 0.15 0.15 0.7 That(j) 10 1.25 1e3];
    for i = 1:length(T)
        opt_R = NaN(length(opt_time),3);
        [t, x] = ode45(@(t,y) mixo_model(t, y, tspan, T(i)*ones(length(tspan),1), I, p, opt_time), tspan, x0, opt);
        u(i,j) = diff(x(end-1:end,8))./x(end-1,8)./diff(t(end-1:end));
    end
end

% Increase parameters
uplus = NaN(length(T),length(That),length(p));
for k = 1:length(p)
    p = params;
    p(13) = 0; % turn off evolution
    p(k) = 1.2*p(k);
    for j = 1:length(That)
        x0 = [1e6 0.15 0.15 0.7 That(j) 10 1.25 1e3];
        for i = 1:length(T)
            opt_R = NaN(length(opt_time),3);
            [t, x] = ode45(@(t,y) mixo_model(t, y, tspan, T(i)*ones(length(tspan),1), I, p, opt_time), tspan, x0, opt);
            uplus(i,j,k) = diff(x(end-1:end,8))./x(end-1,8)./diff(t(end-1:end));
        end
    end
end

% Decrease parameters
uminus = NaN(length(T),length(That),length(p));
for k = 1:length(p)
    p = params;
    p(13) = 0; % turn off evolution
    p(k) = 0.8*p(k);
    for j = 1:length(That)
        x0 = [1e6 0.15 0.15 0.7 That(j) 10 1.25 1e3];
        for i = 1:length(T)
            opt_R = NaN(length(opt_time),3);
            [t, x] = ode45(@(t,y) mixo_model(t, y, tspan, T(i)*ones(length(tspan),1), I, p, opt_time), tspan, x0, opt);
            uminus(i,j,k) = diff(x(end-1:end,8))./x(end-1,8)./diff(t(end-1:end));
        end
    end
end

save(sd6, 'T', 'u', 'uplus', 'uminus');


%% Functions

function dn = mod_diag(x,t,p,T,I)

% QC flag (check for convergence)
dn.qc = var(diff(x(end-10:end,8))./x(end-10:end-1,8)./diff(t(end-10:end))) < 1e-4;

% net growth rate
dn.ngr = diff(x(end-1:end,8))./x(end-1,8)./diff(t(end-1:end));

% strategy
dn.strat = x(end,2:4);

% cellular C
dn.cellc = x(end,6);

% cellular N
dn.celln = x(end,7);

% photosynthesis
dn.photo =  x(end,3)*p(6)*p(7)^((T-p(22))/10)*I/(p(19) + I)/x(end,6);

% grazing rate
dn.graz = x(end,4)*p(8)*p(9)^((T-p(22))/10)*x(end,1)/(p(20) + x(end,1))*p(17)/x(end,6);

% respiration rate
dn.resp = p(10)*(T - x(end,5))^2 + p(11)*p(12)^((T-p(22))/10);

% carbon use efficiency
dn.cue = dn.ngr/(dn.ngr + dn.resp);

end

function dn = ts_diag(x,t,p,T,I)

% time (remove dilution time step)
dn.t = t(~diff(t)==0);

% abundance
dn.x = x(1:end-1,8);
dn.x(diff(t)==0,:) = [];

% thermal optimum
dn.tevo = x(1:end-1,5);
dn.tevo(diff(t)==0) = [];

% net growth rate
dn.ngr = diff(x(:,8))./x(1:end-1,8)./diff(t);
dn.ngr(diff(t)==0) = [];

% strategy
dn.strat = x(1:end-1,2:4);
dn.strat(diff(t)==0,:) = [];

% cellular C
dn.cellc = x(1:end-1,6);
dn.cellc(diff(t)==0) = [];

% cellular N
dn.celln = x(1:end-1,7);
dn.celln(diff(t)==0) = [];

% photosynthesis
dn.photo =  x(1:end-1,3).*p(6).*p(7).^((T(1:end-1)-p(22))./10).*I(1:end-1)./(p(19) + I(1:end-1))./x(1:end-1,6);
dn.photo(diff(t)==0) = [];

% grazing rate
dn.graz = x(1:end-1,4).*p(8).*p(9).^((T(1:end-1)-p(22))./10).*x(1:end-1,1)./(p(20) + x(1:end-1,1)).*p(17)./x(1:end-1,6);
dn.graz(diff(t)==0) = [];

% respiration rate
dn.resp = p(10)*(T(1:end-1) - x(1:end-1,5)).^2 + p(11).*p(12).^((T(1:end-1)-p(22))./10);
dn.resp(diff(t)==0) = [];

% carbon use efficiency
dn.cue = dn.ngr./(dn.ngr + dn.resp);

end
