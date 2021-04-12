clc,clear
close all
tspan = [0 20];
x0 = [13.86 0 0.1 0.1]; % S A X V
[T,C] = ode15s(@batchmodel, tspan, x0);


%%
% parametros
Yas = 0.667; % g/g stochiometric cte
Yoa = 1.067;
Yos = 1.067;
Yxa = 0.4;
Yxsof = 0.15;
Yxsox = 0.51;
Ca = 1/30; % mol Carbon / g acetato
Cs = 1/30; % mol Carbon /g azucar
Cx = 0.04; % mol Carbon / g biomasa
Ka = 0.05; % g/L
KiO = 4; % g/L fitting
KiS = 5; % g/L fitting
Ks = 0.05; % g/L
qAcmax = 0.06; % g/g hr
qm = 0.04; % g/ g hr
qOmax = 15.6; % mmol/g hr
qSmax = 1.3;



%%
qS = ((qSmax)./(1+(C(:,2)./KiS))).*(C(:,1)./(C(:,1)+Ks)); % glucose uptake [g g-1 celulas h-1]
%disp(qS)
qSox = qS;
qSox_an = (qSox - qm).*Yxsox.*Cx/Cs; 
qSox_en = qSox - qSox_an;
qOs = qSox_en.*Yos;
qSof = qS-qSox;
%qSof = qS*0.25;
%overflow
qSof_an = qSof.*Yxsof.*Cx./Cs;
qSof_en = qSof - qSof_an;
qAp = qSof_en.*Yas;
qAc = qAcmax.*C(:,2)./(C(:,2) + Ka);
disp(qAc)
qAc_an = qAc.*Yxa.*Cx./Ca;
qAc_en = qAc - qAc_an;
qO = qOs + qAc_en.*Yoa;
mu = (qSox - qm).*Yxsox + qSof.*Yxsof + qAc.*Yxa;
OCR = (qO/32).*C(:,3).*1000;

%%
figure(1)
subplot(3,2,1)
yyaxis left
plot(T,C(:,1),'g',"LineWidth",2);
ylim([0 15]);
ylabel('Glucose [g/L]')
yyaxis right
plot(T,C(:,3),'r',"LineWidth",2);
ylabel('Biomass [g/L]');
legend('Glucose','Biomass');
ylim([0 10]);
xlabel('Time [h]')

subplot(3,2,2)
plot(T,C(:,2),'k','LineWidth',2);
legend('Acetate')
xlabel('Time [h]')
ylabel('Acetate [g/L]');

subplot(3,2,3)
plot(T,mu,'k','LineWidth',2)
legend('mu')
ylabel('mu [h-1]')
xlabel('Time [h]')

subplot(3,2,4)
plot(T,qAc,'m','LineWidth',2)
legend('qAc')
ylabel('qAc [g g-1 h-1]')
xlabel('Time [h]')

subplot(3,2,5)
plot(T,qO,'k','LineWidth',2)
legend('qO')
xlabel('Time [h]')
ylabel('qO [mmol g-1 h-1]')

subplot(3,2,6)
plot(T,OCR,'g','LineWidth',2)
legend('OCR')
ylabel(' OCR [mmol L-1 h-1]')
xlabel( 'Time [h]')
clear
%%

clc,clear
close all
tspan = [0 12];
x0 = [7.67 0 0.05 0.1]; % S A X V
[T,C] = ode15s(@fedbatch, tspan, x0);

% parma
G = 1.1;
Yas = 0.667; % g/g stochiometric cte
Yoa = 1.067;
Yos = 1.067;
Yxa = 0.4;
Yxsof = 0.15;
Yxsox = 0.51;
Ca = 1/30; % mol Carbon / g acetato
Cs = 1/30; % mol Carbon /g azucar
Cx = 0.04; % mol Carbon / g biomasa
Ka = 0.05; % g/L
KiO = 4; % g/L fitting
KiS = 5; % g/L fitting
Ks = 0.05; % g/L
qAcmax = 0.2; % g/g hr
qm = 0.04; % g/ g hr
qOmax = 13.4; % mmol/g hr
qSmax = 1.25; % g/g hr

qS = ((qSmax)./(1+(C(:,2)./KiS))).*(C(:,1)./(C(:,1)+Ks)); % glucose uptake [g g-1 celulas h-1]
%disp(qS)
qSox = qS;
qSox_an = (qSox - qm).*Yxsox.*Cx/Cs; 
qSox_en = qSox - qSox_an;
qOs = qSox_en.*Yos;
qSof = qS-qSox;
%qSof = qS*0.25;
%overflow
qSof_an = qSof.*Yxsof.*Cx./Cs;
qSof_en = qSof - qSof_an;
qAp = qSof_en.*Yas;
qAc = qAcmax.*C(:,2)./(C(:,2) + Ka);
disp(qAc)
qAc_an = qAc.*Yxa.*Cx./Ca;
qAc_en = qAc - qAc_an;
qO = qOs + qAc_en.*Yoa;
mu = (qSox - qm).*Yxsox + qSof.*Yxsof + qAc.*Yxa;
OCR = (qO/32).*C(:,3).*1000;

% plots

figure(2)
subplot(3,2,1)
yyaxis left
plot(T,C(:,1),'g',"LineWidth",2);
ylim([0 15]);
ylabel('Glucose [g/L]')
yyaxis right
plot(T,C(:,3),'r',"LineWidth",2);
ylabel('Biomass [g/L]');
legend('Glucose','Biomass');
ylim([0 10]);
xlabel('Time [h]')

subplot(3,2,2)
plot(T,C(:,2),'k','LineWidth',2);
legend('Acetate')
xlabel('Time [h]')
ylabel('Acetate [g/L]');

subplot(3,2,3)
plot(T,mu,'k','LineWidth',2)
legend('mu')
ylabel('mu [h-1]')
xlabel('Time [h]')

subplot(3,2,4)
plot(T,qAc,'m','LineWidth',2)
legend('qAc')
ylabel('qAc [g g-1 h-1]')
xlabel('Time [h]')

subplot(3,2,5)
plot(T,qO,'k','LineWidth',2)
legend('qO')
xlabel('Time [h]')
ylabel('qO [mmol g-1 h-1]')

subplot(3,2,6)
plot(T,OCR,'g','LineWidth',2)
legend('OCR')
ylabel(' OCR [mmol L-1 h-1]')
xlabel( 'Time [h]')

