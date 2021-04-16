clc,clear
close all

global  Yax YS_ox_X YS_of_X Yoa K_A C_X C_A C_S qm qO_max qAc_max qS_max mu_set Ysa Yso K_i_A K_S Xin Vin Sfeed

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
qm = 0.04;
qO_max = 13.4*32/1000; % g/ g h
qAc_max = 0.2;
qS_max = 1.25;

mu_set = 0.3;

Xin = 0.4;
Vin = 6.8;
Sfeed = 550;
tspan = [0 24];
x0 = [0.5 0 Xin Vin]; % S A X V O
[T,C] = ode15s(@fedbatch, tspan, x0);

%

qS = (qS_max.*C(:,1)./(K_S + C(:,1))).*1./(1+(C(:,2)./K_i_A));
qOs = min((qS-(qS-qm).*YS_ox_X.*C_X./C_S).*Yso , qO_max);
qSox = ((qOs./Yso)-(qm.*YS_ox_X.*C_X./C_S))./(1-YS_ox_X.*C_X./C_S);
qSof = qS - qSox;
qAp = (qSof - qSof.*YS_of_X.*C_X./C_S).*Ysa;
qAc = min(qAc_max.*C(:,2)./(C(:,2)+K_A) , (qO_max - qOs).*Yoa./(1-Yax.*C_X./C_A));
mu = (qSox - qm).*YS_ox_X + qSof.*YS_of_X + qAc.*Yax;


% oxidative pathway; % glucose uptake [g g-1 celulas h-1]
% relacion rutas metabolicas


%%
figure(1)
subplot(2,2,1)
yyaxis left
plot(T,C(:,1),'g',"LineWidth",2);
ylim([0 15]);
ylabel('Glucose [g/L]')
yyaxis right
plot(T,C(:,3),'r',"LineWidth",2);
ylabel('Biomass [g/L]');
legend('Glucose','Biomass');
xlabel('Time [h]')

subplot(2,2,2)
plot(T,C(:,2),'k','LineWidth',2);
legend('Acetate')
xlabel('Time [h]')
ylabel('Acetate [g/L]');

subplot(2,2,3)
plot(T,C(:,4),'g','LineWidth',2);
legend('Volumen [L]')
xlabel('Time [h]')
ylabel('Volume [L]');

subplot(2,2,4)
plot(T,mu,'g','LineWidth',2)
legend('mu [h-1]')
xlabel('Time [h]')
ylabel('Specific Growth [h-1]')


