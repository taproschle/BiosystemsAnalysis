clc,clear
close all
tspan = [0 20];
x0 = [13.86 0 0.1 0.1 0.035]; % S A X V O
[T,C] = ode15s(@batchmodel, tspan, x0);


%% yields
Yas = 0.23; %  g A /g S stochiometric cte CARCAMO
%Yas = 0.667; % XU
Yoa = 20.66/1000; % g o2 / g A CARCAMO
%Yoa = 1.067; % XU 

Yos = 401/1000; % g o2 / g S CARCAMO
%Yos = 1.067; % XU
%Yxa = 14.55; % g DCW / gIP CARCAMO
Yxa = 0.4; % XU
Yxsof = 0.7; %CARCAMO
%Yxsof = 0.15; % XU
Yxsox = 0.3; % CARCAMO
%Yxsox = 0.51; %XU


%%
Ca = 1/30; % mol Carbon / g acetato
Cs = 1/30; % mol Carbon /g azucar
Cx = 0.04; % mol Carbon / g biomasa

 %%
 %%
%KiO = 4; % g/L fitting
Ka = 0.5236; % CARCAMO
%Ka = 0.05; %XU
KiS = 0.45; % g A /L fitting CARCAMO
%KiS = 5; % XU
Ks = 8; % g S/L CARCAMO
%Ks = 0.05; % XU
Ko = 0.0045; %g O/L CARCAMO
%Ko = 1/100;
qAcmax = 6; % g/g hr CARCAMO
%qAcmax = 0.06; % XU
qm = 0.04; % g/ g hr
qOmax = 83.5/1000; % g/g hr CARCAMO
%qOmax = 15.6*32/1000; % XU
qSmax = 4.9; % g/g hr CARCAMO
%qSmax = 1.3; % Xu
muc = 0.55;
Osat = 0.035;
kla = 180; % [h-1]

% oxidative pathway; % glucose uptake [g g-1 celulas h-1]
% relacion rutas metabolicas
qS_T = qSmax.*C(:,1)./((Ks + C(:,1)).*(1 + C(:,2)./KiS));
qO_crit = qOmax.*C(:,5)./(Ko + C(:,5));
qS_ox = min (qS_T , qO_crit ./Yos) ; 
qS_of   = max(0,qS_T - qS_ox); 
qAc = min(qAcmax.*C(:,2)./(C(:,2) + Ka),(qO_crit - qS_ox.*Yos)./Yoa) ;
% ecs constitutivas 
qS = qS_ox + qS_of;
qA = qS_of.*Yas - qAc; % qA produccion - qA de consumo
mu = qS_ox.*Yxsox + qS_of.*Yxsof + qAc.*Yxa; 
qO2 = qS_ox.*Yos + qAc.*Yoa; % lo q se consume en la via oxidativa + lo que se consume en la via overflow
%mu = (qS_ox-qm).*Yxsox + qS_of.*Yxsof + qAc.*Yxa; 
OCR = (qO2).*C(:,3);

%%
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
ylim([0 10]);
xlabel('Time [h]')

subplot(3,3,2)
plot(T,C(:,2),'k','LineWidth',2);
legend('Acetate')
xlabel('Time [h]')
ylabel('Acetate [g/L]');

subplot(3,3,3)
plot(T,mu,'k','LineWidth',2)
legend('mu')
ylabel('mu [h-1]')
xlabel('Time [h]')

subplot(3,3,4)
plot(T,qAc,'m','LineWidth',2)
legend('qAc')
ylabel('qAc [g g-1 h-1]')
xlabel('Time [h]')

subplot(3,3,5)
plot(T,qO2,'k','LineWidth',2)
legend('qO')
xlabel('Time [h]')
ylabel('qO [gO2 g-1 DW h-1]')

subplot(3,3,6)
plot(T,OCR,'g','LineWidth',2)
legend('OCR')
ylabel(' OCR [gO2 L-1 h-1]')
xlabel( 'Time [h]')

subplot(3,3,7)
plot(T,C(:,5),'g','LineWidth',2)
legend('Oxygen [g/L]')
xlabel(' Time [h]')
ylabel(' Oxygen [g/L]')
clear