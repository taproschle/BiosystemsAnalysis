% Esto codigo corre el Simulink "Lazo_Abierto" para el crecimiento en alta
% densidad celular de Saccharomyces cerevisiae con overflow a 0.14 [1/h]
clc, clear, close all,format compact, warning off, format shortg 

tic
Ts = 1/60; %Tiempo de samplig 1[min]  = 1/60[h]
t_sim = 40;

% Condiciones iniciales del proceso 
x0 = 0.4; s0 = 0.5; a0 = 0; o0 = 0.004; c0 = 0.675; v0 = 0.3;
C_iniciales = [x0  s0  a0   o0   c0   v0]; 
save('CondicionesIniciales.mat ') % Para ser cargadas en la S-function

% Carga de parametros del modelo
ParametrosModelo
V_Stop = 3; % [L]

% Estrategia de control 1) PI,  2) Antireset PID
Estrategia = 1;  
PI_Kp = 1.87082881181595; PI_Ki = 1163.22932064265; PI_Kd = 0;
PI_Kp_2 = 13.0; PI_Ki_2 = 0.06; PI_Kd_2 = 0;
% % ============== Estrategia Rango dividido =========
% Valores para el ruido 
sigma_DO = 10^-8; sigma_RQ= 10^-4;

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
O2s  = A_1_min/Henry;                
spO2 = 0.3*O2s;

% Instruccion para correr el Simulink
sim('Lazo_Abierto')
%% % ======================== Extraccion de resultados =======================
close all
%tiempo
t = Var_Estado(:,11);

% Variables de estado
X = Var_Estado(:,1); S = Var_Estado(:,2); P = Var_Estado(:,3); O = Var_Estado(:,4); C = Var_Estado(:,5); V = Var_Estado(:,6); Scrit_vec = Var_Estado(:,9); F_Alimentacion = Var_Estado(:,10);

% Ecuaciones constitutivas
mu  = Var_Estado(:,7); muc = Var_Estado(:,8);  klao2 = Var_Algebraicas(:,1); klaco2= Var_Algebraicas(:,2); RQ_vec = Var_Algebraicas(:,3); CER = Var_Algebraicas(:,4); OUR = Var_Algebraicas(:,5);

% Variables en ley de control
N = OXIGENO(:,1); G = OXIGENO(:,2); u = OXIGENO(:,3); e = OXIGENO(:,4); yo2 = OXIGENO(:,5); fo2 = OXIGENO(:,6);

% Indices de control
IAE = Indicadores(:,1); ISU = Indicadores(:,2); 
%%
close all
figure (10)
subplot(1,2,1)
plot(t, X,'k','LineWidth',1.5) ;grid; ylabel('X [gDCW/L]','FontSize',14);xlabel('t [h]','FontSize',14)

subplot(1,2,2)
yyaxis left
plot(t, V,'LineWidth',1.5); grid on; ylabel('V [L]','FontSize',14); xlabel('t [h]','FontSize',14)
yyaxis right
plot(t,F_Alimentacion*1000/60,'LineWidth',1.5); grid on; ylabel('F [mL/min]','FontSize',14); xlabel('t [h]','FontSize',14)

%%
% ===================== Figuras de resultados =====================
figure (1)
subplot(3,3,1)
plot(t, X,'k','LineWidth',1.5) ;grid; ylabel('X [gDCW/L]','FontSize',14);xlabel('t [h]','FontSize',14)

subplot(3,3,2)
yyaxis left
plot(t, V); grid on; ylabel('V [L]','FontSize',14); xlabel('t [h]','FontSize',14)
yyaxis right
plot(t,F_Alimentacion); grid on; ylabel('F [L/h]','FontSize',14); xlabel('t [h]','FontSize',14)

subplot(3,3,3)
plot(t, C,'r','LineWidth',1.5); grid; ylabel('DCO_{2} [g/L]','FontSize',14);xlabel('t [h]','FontSize',14)

subplot(3,3,4)
plot(t, S,'m',t,Scrit_vec,'r','LineWidth',1.5);grid; ylabel('S [g/L]','FontSize',14);xlabel('t [h]','FontSize',14)
legend('S' ,' S_{crit}','Location','best','FontSize',10)

subplot(3,3,5); 
plot(t, mu,'m',t, muc,'r','LineWidth',1.5); grid; ylabel('\mu [1/h]','FontSize',14);xlabel('t [h]','FontSize',14)
legend('\mu' ,' \mu_{crit}','Location','best','FontSize',10)

subplot(3,3,6)
plot(t, P,'r','LineWidth',1.5); grid; ylabel('P [g/L]','FontSize',14); xlabel('t [h]','FontSize',14)

subplot(3,3,7)
plot(t, RQ_vec,'m','LineWidth',1.5); ylabel('RQ [mol CO_2/mol O_2]','FontSize',14); grid on; xlabel('t [h]','FontSize',14)
xlim([0.1, max(t)])

spO2_vec = linspace(spO2,spO2,length(t))*1000;
subplot(3,3,8)
plot(t, O*1000,'g',t, spO2_vec,'k--', 'LineWidth',1.5); grid; ylabel('DO [mg/L]','FontSize',14);xlabel('t [h]','FontSize',14)

subplot(3,3,9)
plot(t, yo2,'g','LineWidth',1.5); ylabel('y_{O2}','FontSize',14); grid; xlabel('t [h]','FontSize',14)

figure(2)
subplot(2,3,1)
plot(t, klao2,t, klaco2,'LineWidth',1.5); grid; ylabel('k_La [1/h]','FontSize',14);xlabel('t [h]','FontSize',14)
legend('k_La_{O2}' ,' k_La_{CO2}','Location','best','FontSize',10)

subplot(2,3,2)
plot(t, N,'g','LineWidth',1.5); ylabel('N [rpm]','FontSize',14); grid; xlabel('t [h]','FontSize',14)

subplot(2,3,3)
plot(t, u*100,'m','LineWidth',1.5); ylabel('u [%]','FontSize',14); grid;xlabel('t [h]','FontSize',14)

subplot(2,3,4)
plot(t,OUR,t, CER,'LineWidth',1.5); ylabel(' [mol/L/h]','FontSize',14); grid; xlabel('t [h]','FontSize',14)
legend('OUR','CER','Location','best','FontSize',10)

subplot(2,3,5)
plot(t, G,'r','LineWidth',1.5); ylabel('G [L/min]','FontSize',14); grid;xlabel('t [h]','FontSize',14)

subplot(2,3,6)
plot(t, e*1000,'k','LineWidth',1.5); grid; ylabel('Error DO [mg/L]','FontSize',14); xlabel('t [h]','FontSize',14)
%%
% ====================== Indices de control ====================

disp('                              Ïndices de control del lazo DO')
disp('indice_IAE')
indice_e2 = sum(IAE);
disp (indice_e2)

disp('indice_ISU')
[k,l,indice_u2] = find(ISU,1,'last');
disp (indice_u2)
toc