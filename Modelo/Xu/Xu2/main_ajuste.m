%% Ajuste parámetros XU
clc; clear; close all;
tabla = readtable('datos_troles.csv');
tabla = table2array(tabla);

texp = tabla(:,1)';
yexp = tabla(:,2:6);

tsim = texp(end);

Yax = 0.667; YS_ox_X = 0.51; YS_of_X = 0.15; Yoa = 1.067; Ysa = 0.667; 
Yso = 1.067; C_X = 0.04; C_A = 1/30; C_S = 1/30; K_A = 0.05; K_i_A = 5;
K_S = 0.05; qm = 0.04; qO_max = 13.4*32/1000; qAc_max = 0.2; qS_max = 1.25;
x0 = [YS_ox_X YS_of_X Yax Yoa Ysa Yso C_X C_A C_S qm ...
    qS_max qAc_max qO_max K_S K_i_A K_A];

% MEIGO Settings
clc
problem.f = 'funObj';
problem.x_L = 0*ones(1,16);
xU = 2*ones(1,16); xU(15) = 10;
problem.x_U = xU;
problem.x_0 = x0;

opts.maxeval = 100;
opts.maxtime = 1000;
opts.iterprint = 1;

opts.ndiverse = 100;
opts.local.solver = 'fminsearch';
opts.local.finish = 'solnp';

id = 'MATLAB:ode15s:IntegrationTolNotMet';
warning('off',id)

Results = MEIGO(problem, opts, 'ESS');

params = Results.xbest;

%% 
k = params;
% Parámetros no ajustables:
mu_set = 0.13;
klao2 = 180*100;
Xin = 5; Sin = 0.04; Ain = 0; Oin = 4e-3; Vin = 0.3;
Sfeed = 550;
O_sat = 0.035; %850/1000;
K_O =  0.0001; % g o2 L-1

v = [mu_set klao2 Vin Xin Sfeed O_sat K_O];

tspan = [0 tsim];
y0 = [5 0.04 0 0.004 0.3];
fun = @(t,y) xu_model(t,y,v,k);
[T,C] = ode23s(fun, tspan, y0);

% Gráfico:
figure(1)
subplot(2,2,1)
yyaxis left
plot(T,C(:,2),'g',"LineWidth",2);
hold on
plot(texp,yexp(:,2),'og','linewidth',1)
hold off
ylabel('Glucose [g/L]')
yyaxis right
plot(T,C(:,1),'r',"LineWidth",2);
hold on
plot(texp,yexp(:,1),'or','linewidth',1)
hold off
ylabel('Biomass [g/L]');
legend('Glucose','Biomass');
xlabel('Time [h]')
xlim([0 tsim])

subplot(2,2,2)
plot(T,C(:,3),'k','LineWidth',2);
hold on
plot(texp,yexp(:,3),'ok','linewidth',1)
hold off
legend('Acetate')
xlabel('Time [h]')
ylabel('Acetate [g/L]');
xlim([0 tsim])

subplot(2,2,3)
plot(T,C(:,5),'g','LineWidth',2);
hold on
plot(texp,yexp(:,5),'og','linewidth',1)
hold off
legend('Volumen [L]')
xlabel('Time [h]')
ylabel('Volume [L]');
xlim([0 tsim])

subplot(2,2,4)
plot(T,C(:,4),'b','LineWidth',2)
hold on
plot(texp,yexp(:,2),'ob','linewidth',1)
hold off
legend('Oxygen [g/L]')
xlabel('Time [h]')
ylabel('Oxygen [g/L]')
xlim([0 tsim])

var = ["YS_ox_X" , "YS_of_X" , "Yax" , "Yoa" "Ysa" , "Yso", "C_X", ...
"C_A", "C_S", "qm","qS_max", "qAc_max", "qO_max", "K_S", "K_i_A", "K_A"];
tab = table(var',x0',k');

%% Análisis de residuos
% no funca
% options = optimset('MaxIter', 1, 'MaxFunEvals', 1);
% 
% [x,~,residuals,~,~,~,J] = lsqcurvefit(@avoidCurveFit, k, texp, yexp, [],[],options);
% 
% [~, CV] = intconfianza(J, residuals, k, 0.05);
% tvalue = 1./CV * 100;
% display(tvalue)
% 
% fun = @(t,y,k) xu_model(t,y,v,k)
% 
% clc
% identifica((1:2),params,[2 30], data.time, fun , 0.95)
% 
% close all
% ksensibilidad((1:2),params,[2 30],data.time,@odeW13,1)





