% Ajuste de parametros
clear ; clc ; close all
tic
load data.csv

texp    = data(:,1)';
yexp    = data(:,2:4);

tsim    = texp(end);
%          qSmax Ysoxx    Yse  Ysofx Kec
x0      = [3.996 0.296 0.9634 0.3294   1 1]; 
                
% MEIGO Settings
problem.f = 'funObj';
x_B = [ 0.1   100;      % qSmax
        0.2    1;      % Ysoxx
        0.1    100;      % Yse
        0.2   1        % Ysofx
        0.1   100      % Kec  
        1   100];    % Kio

    
problem.x_L = x_B(:,1);
problem.x_U = x_B(:,2);
problem.x_0 = x0;

opts.maxeval = 1000;
opts.maxtime = 100;
opts.iterprint =  1;

opts.ndiverse = 2000;
opts.local.solver = 'fminsearch';
opts.local.finish = 'fminsearch';

id = 'MATLAB:ode15s:IntegrationTolNotMet';
warning('off',id)

Results = MEIGO(problem, opts, 'ESS');

params  = Results.xbest;

%%

k = params;
save('Results\kDew.mat','k');

% Fixed parameters
muset   = 0.197;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 2.4e-3;
Ko      = 0.0001;
v       = [muset X0 V0 Sin klao2 osat Ko];

% Initial conditions
S0 = 0.001;
E0 = 3.95;
O0 = osat;
y0 = [X0 S0 E0 O0 V0];

tspan   = [0 tsim];
fun = @(t,y) dewasme_unified(t,y,v,k);
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'NonNegative',[1,2,3,4,5]);
[T,C] = ode15s(fun,tspan,y0,options);
tDew = T; cDew = C;
Dewdata = [tDew cDew];
save('Results\Dewdata.mat','Dewdata')

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
ylabel('Volume (L)')
xlabel('Time (h)')

toc

