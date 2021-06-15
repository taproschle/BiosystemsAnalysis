% Ajuste de parametros
clear ; clc ; close all
load data.csv

texp    = data(:,1)';
yexp    = data(:,2:5);

tsim    = texp(end);

x0      = [ 0.0500   % Ks
            1.2500   % qSmax
            0.5100   % Ysoxx
            0.0040   % qm
            1.0670   % Yos
            5.0000   % Kie
            0.6670   % Yes
            0.0500   % Kec
            0.2000   % qEmax
            0.1500   % Ysofx
            1.0670   % Yoe
            0.4000   % Yxe
            0.4300]; % qOmax

% MEIGO Settings
problem.f = 'funObj';
x_L = [1e-2 1e-2 0.1 1e-4 1e-2 1 1e-2 1e-3 1e-2 1e-3 1e-3 1e-3 1e-3];
x_U = [1 2 0.8 0.1 2 10 2 1 2 1 2 1 0.4];
problem.x_L = x_L;
problem.x_U = x_U;
problem.x_0 = x0;

opts.maxeval = 1000;
opts.maxtime = 100;
opts.iterprint = 1;

opts.ndiverse = 2000;
opts.local.solver = 'fmincon';
opts.local.finish = 'solnp';

id = 'MATLAB:ode15s:IntegrationTolNotMet';
warning('off',id)

Results = MEIGO(problem, opts, 'ESS');

params  = Results.xbest;

%%

k = params;
save('kXu.mat','k');



% Initial conditions




% Fixed parameters
muset   = 0.11;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.008;
Ko      = 0.0001;
Kio     = 4;
v       = [muset X0 V0 Sin klao2 osat Ko Kio];
save('vXu.mat','v');

% Adjusted parameters (overflow)
% Kie     = 5;
% Yes     = 0.667;
% Kec     = 0.05;
% qEmax   = 0.2;
% Ysofx   = 0.15;
% Yoe     = 1.067;
% Yxe     = 0.4;
% qOmax   = 0.43;
% Kio     = 4;
% kof     = [Kie Yes Kec qEmax Ysofx Yoe Yxe qOmax Kio];

% Initial conditions
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;
y0 = [X0 S0 E0 O0 V0];

tspan   = [0 tsim];
fun = @(t,y) xu_unified(t,y,v,k);
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'NonNegative',[1,2,3,4,5]);
[T,C] = ode15s(fun,tspan,y0,options);
tXu = T; cXu = C;
Xudata = [tXu cXu];
save('Xudata.mat','Xudata')

c1 = "#1B9E77"; c2 = "#D95F02"; c3 = "#7570B3"; c4 = "#E7298A";
c5 = "#66A61E"; c6 = "#E6AB02"; c7 = "#A6761D"; c8 = "#666666";

tiledlayout(2,3)

nexttile
plot(T,C(:,1),'Color',c1,'LineWidth',1.5)
grid on
hold on
plot(texp,yexp(:,1),'s','Color',c1,'LineWidth',1)
ylabel('Biomass (g/L)')
xlabel('Time (h)')

nexttile
plot(T,C(:,2),'Color',c2,'LineWidth',1.5)
grid on
hold on
plot(texp,yexp(:,2),'s','Color',c2,'LineWidth',1)
ylabel('Glucose (g/L)')
xlabel('Time (h)')

nexttile
plot(T,C(:,3),'Color',c3,'LineWidth',1.5)
grid on
hold on
plot(texp,yexp(:,3),'s','Color',c3,'LineWidth',1)
ylabel('Ethanol (g/L)')
xlabel('Time (h)')

nexttile
plot(T,C(:,4),'Color',c4,'LineWidth',1.5)
grid on
hold on
plot(texp,yexp(:,4),'s','Color',c4,'LineWidth',1)
ylabel('Dissolved Oxygen (g/L)')
xlabel('Time (h)')

nexttile
plot(T,C(:,5),'Color',c5,'LineWidth',1.5)
grid on
hold on
% plot(texp,yexp(:,5),'s','Color',c5,'LineWidth',1)
ylabel('Volume (L)')
xlabel('Time (h)')


