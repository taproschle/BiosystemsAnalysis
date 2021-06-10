% Ajuste de parametros
clear ; clc ; close all
load data.csv

texp    = data(:,1)';
yexp    = data(:,2:5);

tsim    = texp(end);

x0      = [ 0.0370   % Ks
            0.6356   % qSmax
            0.5333   % Ysoxx
            0.0129   % qm
            1.5620   % Yos
            1.2399   % Kie
            0.2268   % pEmax
            0.5052   % Kep
            0.9097   % Yes
            0.0134   % Kec
            0.1148   % qEmax
            2.1231   % Kis
            0.2268   % Ysofx
            0.5440   % Yoe
            0.5718]; % Yxe
        
% MEIGO Settings
problem.f = 'funObj';
problem.x_L = 1e-6*ones(1,15);
problem.x_U = 20*ones(1,15);
problem.x_0 = x0;

opts.maxeval = 1000;
opts.maxtime = 1000;
opts.iterprint = 1;

opts.ndiverse = 2000;
opts.local.solver = 'lsqnonlin';
opts.local.finish = 'solnp';

id = 'MATLAB:ode15s:IntegrationTolNotMet';
warning('off',id)

Results = MEIGO(problem, opts, 'ESS');

params  = Results.xbest;

%%

k = params;
save('kAn.mat','k');

% Fixed parameters
muset   = 0.11;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.035;
Ko      = 0.0001;
v       = [muset X0 V0 Sin klao2 osat Ko];

% % Adjusted parameters (overflow)
% Kie     = 1.2399;
% pEmax   = 0.2268;
% Kep     = 0.5052;
% Yes     = 0.9097;
% Kec     = 0.0134;
% qEmax   = 0.1148;
% Kis     = 2.1231;
% Ysofx   = 0.2268;
% Yoe     = 0.5440;
% Yxe     = 0.5718;
% kof     = [Kie pEmax Kep Yes Kec qEmax Kis Ysofx Yoe Yxe];

% Initial conditions
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;
y0 = [X0 S0 E0 O0 V0];

tspan   = [0 tsim];
fun = @(t,y) anane_unified(t,y,v,k);
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'NonNegative',[1,2,3,4,5]);
[T,C] = ode15s(fun,tspan,y0,options);
tAn = T; cAn = C;
Andata = [tAn cAn];
save('Andata.mat','Andata')

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


