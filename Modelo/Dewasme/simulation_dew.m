%% Dewasme Simulation
clc; clear; close all;

tsim = 40;
Ts = 1/60; %Tiempo de samplig 1[min]  = 1/60[h]

% Estrategia de control 1) PI,  2) Antireset PID
Estrategia = 1;  
PI_Kp = 1.87082881181595; PI_Ki = 1163.22932064265; PI_Kd = 0;
PI_Kp_2 = 13.0; PI_Ki_2 = 0.06; PI_Kd_2 = 0;

% % ============== Estrategia Rango dividido =========
% Valores para el ruido 
sigma_DO =  10^-8; % Activar o desactivar ruido en el sensor
sigma_RQ=   10^-4;

% ===================== Valores de upper y lower de N y G y control =======================

% Saturador principal
u_1_low = 0;  u_1_upp = 1;

% Saturador Agitacion [% y rpm]
N_1_low = 0; N_1_upp = 0.30;
N_1_min = 300;   N_1_max = 300; % N total es la suma del min y el max

% Saturador flujo de aire [% y LPM]
% (u(1)-G_1_low)*G_1_mid/G_1_min+G_1_min
G_1_low = 0.30; G_1_upp = 0.50;   
G_1_min = 0.2; G_1_max = 0.6;
G_1_mid = G_1_max - G_1_min;

% Saturador fraccion volumetrica de oxigeno yO2 [% y %v/v]
% (u(1)-A_1_low)*A_1_mid/A_1_low+A_1_min
A_1_low = 0.50; A_1_upp = 1;

% Fraciones volumetricas de oxigeno
A_1_min = 0.21;   % fraccion de oxigeno en aire
A_1_max = 1;    % fraccion máxima deseada - Sin perturbación hasta 0.62
A_1_mid = A_1_max - A_1_min;

% ==================  Valor de KlaO2 (Filtro y Bloque de Proceso) ====
Henry = 26.409;
O2s  = A_1_min/Henry;                
spO2 = 0.3*O2s;

% PARÁMETROS:
% load kDew.mat k
load vDew.mat v

% Ajustables
% Ks      = k(1);
% qSmax   = k(2);
% Ysoxx   = k(3);
% Yso     = k(4);
% Kio     = k(5);
% Yse     = k(6);
% Kec     = k(7);
% Ysofx   = k(8);
% Yeo     = k(9);
% Yex     = k(10);
% qOmax   = k(11);
% Yosof   = k(12);

% Fijos (Cambiar dependiendo el ajuste)
mu_set  = 0.16;
Xin     = 4.21;
Vin     = 1.5;
Sfeed   = 450;
Ko      = 0.0001;

% Seleccionar que datos cargar
% load data_ox.csv
load data_of.csv
t = data_of(:,1);
X = data_of(:,2);
S = data_of(:,3);
E = data_of(:,4);
O = data_of(:,5);

% Correr hasta aquí para hacer el análisis de sensibilidad y ajuste en
%  simulink

% En base a lo realizado por la sesión de ajuste en simulink
% Ajuste con regimen oxidativo
Ks      = 0.155434749;
qSmax   = 4.39560768;
Ysoxx   = 0.325626329;
Yso     = 0.309409692;
Kio     = 5.2044717434937;
Yse     = 1.030307418;
Kec     = 0.005523424;
Ysofx   = 0.317896552;
Yeo     = 2.180279923;
Yex     = 1.928401808;
qOmax   = 0.384622039;
Yosof   = 0.00011;

%% SIMULACION

simulation = sim('dew_sim');

T = simulation.time;
data = simulation.data;
inlets = simulation.inlets;
Indicadores = simulation.Indicadores;
IAE = Indicadores(:,1); ISU = Indicadores(:,2);
F = inlets(:,1); N = inlets(:,2); G = inlets(:,3); yO2 = inlets(:,4);
X = data(:,1); S = data(:,2); A = data(:,3); O = data(:,4);
V = data(:,5); mu = data(:,6); mucrit = data(:,7); Scrit = data(:,8);

dataDew = [T data inlets];

save('simdata_dew','dataDew')

%% Gráfico:

load data_ox.csv

texp    = data_ox(:,1)';
yexp    = data_ox(:,2:5);

figure(1)
subplot(2,2,1)
yyaxis left
plot(T,S,'g',"LineWidth",2); hold on;
plot(texp,yexp(:,2),'og',"LineWidth",1); hold off;
ylabel('Glucose [g/L]')
yyaxis right
plot(T,X,'r',"LineWidth",2); hold on;
plot(texp,yexp(:,1),'or',"LineWidth",1); hold off;
ylabel('Biomass [g/L]');
legend('Glucose','Exp. Glucose','Biomass','Exp. Biomass','location','northwest');
xlabel('Time [h]')
xlim([0 tsim])

subplot(2,2,2)
plot(T,A,'k','LineWidth',2); hold on;
plot(texp,yexp(:,3),'ok','linewidth',1); hold off
legend('Ethanol','Exp. Ethanol')
xlabel('Time [h]')
ylabel('Ethanol [g/L]');
xlim([0 tsim])

subplot(2,2,3)
plot(T,V,'g','LineWidth',2);
legend('Volumen [L]','location','northwest')
xlabel('Time [h]')
ylabel('Volume [L]');
xlim([0 tsim])

subplot(2,2,4)
plot(T,O,'g','LineWidth',2); hold on;
plot(texp,yexp(:,4),'og','linewidth',1); hold off
legend('Oxygen','Exp. Oxygen')
xlabel('Time [h]')
ylabel('Oxygen [g/L]')
xlim([0 tsim])
%% Growth plots
figure(2);
plot(T,mu,'-k','linewidth',2); hold on;
plot(T,mucrit,'-r','linewidth',2); hold off;
legend('\mu','\mu_{crit}'); grid on;
xlabel('Time [h]')
ylabel('Growth Constant [g/L]')


%% Inlet Plots
figure(3);
subplot(2,2,1)
plot(T,F,'-k','linewidth',2)
xlim([0 tsim])
ylabel('F [l/hr]')
xlabel('t [hr]')
legend('Inlet Feed')

subplot(2,2,2)
plot(T,N,'-r','linewidth',2)
xlim([0 tsim])
ylabel('N [rpm]')
xlabel('t [hr]')
legend('Agitation rate')

subplot(2,2,3)
plot(T,G,'-g','linewidth',2)
xlim([0 tsim])
ylabel('G [l/min]')
xlabel('t [hr]')
legend('Inlet Air flow','location','southeast')

subplot(2,2,4)
plot(T,100*yO2,'-b','linewidth',2)
xlim([0 tsim])
ylabel('yO2')
xlabel('t [hr]')
legend('Inlet O_2 %','location','northwest')


%%
% ====================== Indices de control ====================

disp('========== Ïndices de control del lazo DO ===========')
disp('indice_IAE')
indice_e2 = sum(IAE);
disp (indice_e2)

disp('indice_ISU')
[k,l,indice_u2] = find(ISU,1,'last');
disp (indice_u2)
disp('======================================================')








