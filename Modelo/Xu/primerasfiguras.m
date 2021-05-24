clc,clear
close all

global O_sat klao2 K_O Yax YS_ox_X YS_of_X Yoa K_A C_X C_A C_S qm qO_max qAc_max qS_max mu_set Ysa Yso K_i_A K_S Xin Vin Sfeed

O_sat = 0.035;
klao2 = 180*10;
Yax = 0.667;
YS_ox_X = 0.51;
YS_of_X = 0.15;
Yoa = 1.067;
Ysa = 0.667; 
Yso = 1.067;
C_X = 0.04;
C_A = 1/30;
C_S = 1/30;
K_A = 0.05;
K_i_A = 5;
K_S = 0.05;
K_O =  0.0001; % g o2 L-1
qm = 0.04;
qO_max = 13.4*32/1000; % g/ g h
qAc_max = 0.2;
qS_max = 1.25;

mu_set = 0.3;

Xin = 0.4;
Vin = 6.8;
Sfeed = 550;
tspan = [0 24];
x0 = [0 0 Xin Vin 0.035]; % S A X V O
%x0 = [0 0 Xin Vin]; % S A X V O
[T,C] = ode15s(@fedbatch, tspan, x0);

%

qS = (qS_max.*C(:,1)./(K_S + C(:,1))).*1./(1+(C(:,2)./K_i_A));
%qOs = min((qS-(qS-qm)*YS_ox_X*C_X/C_S)*Yso , qO_max); 

%qOs = min(qO_max*(O/(K_O+O))*(1/(1+(A/K_i_A))),qO_max*(1/(1+(A/K_i_A))));
qOs = min(qO_max.*(C(:,4)./(K_O+C(:,4))).*(1./(1+(C(:,2)./K_i_A))),qO_max);
qSox = min(((qOs/Yso)-(qm*YS_ox_X*C_X/C_S))/(1-YS_ox_X*C_X/C_S),qS);
%qSox = min(qOs./Yso,qS);
%qSox = qOs/Yso;
qSof = max(qS - qSox,0);

qAp = (qSof - qSof.*YS_of_X.*C_X./C_S).*Ysa;
qAc = min(qAc_max.*C(:,2)./(C(:,2)+K_A) , (qO_max - qOs).*Yoa./(1-Yax.*C_X./C_A));
mu = (qSox - qm).*YS_ox_X + qSof.*YS_of_X + qAc.*Yax;
qO = qOs + (qAc - qAc.*Yax.*C_X./C_A).*Yoa;
    

% oxidative pathway; % glucose uptake [g g-1 celulas h-1]
% relacion rutas metabolicas


% %
figure(1)
subplot(3,3,1)
yyaxis left
plot(T,C(:,1),'g',"LineWidth",2);
ylim([0 15]);
ylabel('Glucose [g/L]')
yyaxis right
plot(T,C(:,3),'r',"LineWidth",2);
ylabel('Biomass [g/L]');
legend('Glucose','Biomass');
xlabel('Time [h]')

subplot(3,3,2)
plot(T,C(:,2),'k','LineWidth',2);
legend('Acetate')
xlabel('Time [h]')
ylabel('Acetate [g/L]');

subplot(3,3,3)
plot(T,C(:,4),'g','LineWidth',2);
legend('Volumen [L]')
xlabel('Time [h]')
ylabel('Volume [L]');

subplot(3,3,4)
plot(T,mu,'g','LineWidth',2)
legend('mu [h-1]')
xlabel('Time [h]')
ylabel('Specific Growth [h-1]')

subplot(3,3,5)
plot(T,C(:,5),'g','LineWidth',2)
legend('Oxygen [g/L]')
xlabel('Time [h]')
ylabel('Oxygen [g/L]')

subplot(3,3,6)
plot(T,qO,'g','LineWidth',2)
legend('oxygen consumption rate [h-1]')
xlabel('Time [h]')
ylabel('qO')

subplot(3,3,7)
plot(T,qS,'g','LineWidth',2)
legend('glucose consumption rate [h-1]')
xlabel('Time [h]')
ylabel('qS')

subplot(3,3,8)
plot(T,qSox,'g','LineWidth',2)
legend('glucose consumption rate in OP [h-1]')
xlabel('Time [h]')
ylabel('qSox')

