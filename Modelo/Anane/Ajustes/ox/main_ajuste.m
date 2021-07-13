% Ajuste de parametros
clear ; clc ; close all
load data.csv

texp    = data(:,1)';
yexp    = data(:,2:5);

tsim    = texp(end);

  
x0 = [0.6356   % qSmax
      0.1148   % qEcmax
      0.4267   % Kec
      0.5333   % Ysoxx
      0.2268   % Ysofx
      0.5718]; % Yxe
        
% MEIGO Settings
problem.f = 'funObj';
x_B =  [1e-5, 100;   % qSmax
        1e-5, 100;   % qEmax
        1e-5, 1;     % Kec     
        1e-5, 0.8;   % Ysoxx
        1e-5, 1;     % Ysofx
        1e-5, 1];    % Yxe
    
problem.x_L = x_B(:,1);
problem.x_U = x_B(:,2);
problem.x_0 = x0;

opts.maxeval = 1000;
opts.maxtime = 1000;
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
save('Results\kAn.mat','k');

% Fixed parameters
muset   = 0.11;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.008;
Ko      = 0.0001;
v       = [muset X0 V0 Sin klao2 osat Ko];

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
% save('Andata.mat','Andata')

c1 = "#1B9E77"; c2 = "#D95F02"; c3 = "#7570B3"; c4 = "#E7298A";
c5 = "#66A61E"; c6 = "#E6AB02"; c7 = "#A6761D"; c8 = "#666666";

figure(1)
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
ylabel('Volume (L)')
xlabel('Time (h)')


