%% Xu Simulation
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
% load kAn.mat k
% load vAn.mat v
% mu_set  = v(1);
% Xin     = v(2);
% Vin     = v(3);
% Sfeed   = v(4);
% klao2   = v(5);
% osat    = v(6);
% Ko      = v(7);

mu_set   = 0.11;
Xin      = 4.125;
Vin      = 0.3;
Sfeed     = 450;
klao2   = 180*100;
osat    = 0.008;
Ko      = 0.0001;
v       = [mu_set Xin Vin Sfeed klao2 osat Ko];
save('vAn.mat','v')

% Parametros ajustados por simulink (Fijando por FIA)
Kep = 0.73784;
Kie = 2.4227;
Kis = 1.8205;
Ks = 0.0044503;
Yeo = 0.93296;
Yse = 0.10329;
Yso = 0.8125;
qEpmax = 0.73184;
qm = 0.0068952;
Kec = 0.29725;
Yex = 0.50463;
Ysofx = 0.23879;
Ysoxx = 0.2886;
qEcmax = 0.32776;
qSmax = 2.0811;

% Seleccionar que datos cargar
data = table2array(readtable('data_ox.csv'));
t = data(:,1);
X = data(:,2);
S = data(:,3);
E = data(:,4);
O = data(:,5);


%%
% SIMULACION

simulation = sim('an_sim');

T = simulation.time;
Data = simulation.data;
inlets = simulation.inlets;
Indicadores = simulation.Indicadores;
IAE = Indicadores(:,1); ISU = Indicadores(:,2);
F = inlets(:,1); N = inlets(:,2); G = inlets(:,3); yO2 = inlets(:,4);
X = Data(:,1); S = Data(:,2); A = Data(:,3); O = Data(:,4); 
V = Data(:,5);

dataAn = [T Data inlets];

save('simdata_an','dataAn')

%% Gráfico:

load data_ox.csv

texp    = data_ox(:,1)';
yexp    = data_ox(:,2:4);

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
legend('Acetate','Exp. Acetate')
xlabel('Time [h]')
ylabel('Acetate [g/L]');
xlim([0 tsim])

subplot(2,2,3)
plot(T,V,'g','LineWidth',2);
legend('Volumen [L]','location','northwest')
xlabel('Time [h]')
ylabel('Volume [L]');
xlim([0 tsim])

subplot(2,2,4)
plot(T,O,'g','LineWidth',2); hold on;
legend('Oxygen')
xlabel('Time [h]')
ylabel('Oxygen [g/L]')
xlim([0 tsim])

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
[K,l,indice_u2] = find(ISU,1,'last');
disp (indice_u2)
disp('======================================================')








