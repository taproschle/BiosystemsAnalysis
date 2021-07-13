% Ajuste de parametros
clear ; clc ; close all

%load dataferm2.csv
load data_ox.csv;

texp    = data_ox(:,1)';
yexp    = data_ox(:,2:4);

tsim    = texp(end);
%               qSmax     qOmax     Ysoxx    Ysofx   Yso
x0      = [  2      0.1            0.1         0.5             0.5 ];  
 % PARAMETROS  qsmax 42.6082 qomax 0.239492 Ysoxx 0.2118  Ysofx 0.730404
 % Yso 0.846497


% MEIGO Settings
problem.f = 'funObj';
x_L = [ 1e-3      1e-3         1e-2     	1e-4        1e-2];
x_U = [ 100    10         1          1           5];
problem.x_L = x_L;
problem.x_U = x_U;
problem.x_0 = x0;

opts.maxeval = 1000;
opts.maxtime = 100;
opts.iterprint = 1;

opts.ndiverse = 2000;
opts.local.solver = 'fminsearch';
opts.local.finish = 'fminsearch';

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
X0      = 4.15;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.008;
Ko      = 0.0001;
Kio     = 4;
v       = [muset X0 V0 Sin klao2 osat Ko Kio];
save('vXu.mat','v');


% Initial conditions
S0 = 0.001;
E0 = 4.1;
O0 = osat*0.3;
y0 = [X0 S0 E0 O0 V0];
% y0 = [X0 S0 E0 V0];


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
ylabel('Dissolved Oxygen (g/L)')
xlabel('Time (h)')

nexttile
plot(T,C(:,5),'Color',c5,'LineWidth',1.5)
grid on
hold on
% plot(texp,yexp(:,5),'s','Color',c5,'LineWidth',1)
ylabel('Volume (L)')
xlabel('Time (h)')


