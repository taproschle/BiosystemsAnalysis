clear,clc
close all
global Vin Xin Factor_mu YS_ox_X YS_ox_O YS_f_X YS_f_IP YIP_X YIP_O qS_max qIP_max qO_max K_S K_IP Ki_IP K_O kLa_O2 Sfeed DO_max mu_max mu_set 



YS_ox_X = 0.3;
YS_ox_O = 401;
YS_f_X = 0.7;
YS_f_IP = 0.23;
YIP_X = 14.55;
YIP_O = 20.66;
qS_max = 4.9;
qIP_max = 6;
qO_max = 83.5; 
K_S = 8;
K_IP = 0.45;
Ki_IP = 10;
K_O = 4.5;
kLa_O2 = 140*20;


Vin =  6.8;
Xin = 0.4;
Sfeed = 500;
Factor_mu =0;
mu_set = 0.3;
mu_max = 0.65;
% mu_set = 0.2;
% u = mu_set; % mu_Set XIN Vin Sin
tspan = [0 20];

x0 = [Xin 0.1 0 Vin 35]; % X S Et V DO

[t1,x1] = ode15s(@(t,x) modelo_carcamo(t,x), tspan, x0); % X S P V DO

 qS_crit = qO_max ./ YS_ox_O ; %  [g_Glu / (g_DCW * h)]
% tasa específica de crecimiento crítica
% mu_crit = qS_crit .* YS_ox_X ; % [1/h]
% mu_set = mu_crit .* Factor_mu;



F = mu_set.*Vin.*Xin.*exp(mu_set.*t1)./(YS_ox_X.*Sfeed);



q_ST = qS_max.*x1(:,2) ./ ((K_S + x1(:,2)).*(1 + x1(:,3)./Ki_IP));
qO_crit = qO_max*x1(:,5) ./ (K_O + x1(:,5));

qS_ox = min (q_ST , qO_crit ./ YS_ox_O ) ; 
qS_f   = max (0 , q_ST - qS_ox) ; 
qIP_cons = min (qIP_max .* x1(:,3) ./ (x1(:,3) + K_IP) , (qO_crit - qS_ox.*YS_ox_O) ./ YIP_O) ;

% Ecuaciones constitutivas 

qS = qS_ox + qS_f;
qIP = qS_f.*YS_f_IP - qIP_cons;
mu = min(mu_max,qS_ox.*YS_ox_X + qS_f.*YS_f_X + qIP_cons.*YIP_X); 
%mu = qS_ox.*YS_ox_X + qS_f.*YS_f_X + qIP_cons.*YIP_X;
qO2 = qS_ox.*YS_ox_O + qIP_cons.*YIP_O;

figure(1)
subplot(2,2,1)
plot(t1,x1(:,1),'g',t1,x1(:,2),'r',t1,x1(:,3),'k','LineWidth',2);
legend('Biomass','Substract','Ethanol')
xlabel(' Time (h)')
ylabel(' Concentration [g/L] ');



subplot(2,2,2);
plot(t1,F,'r')
ylabel('Flujo')


subplot(2,2,3)
plot(t1,mu,'r')

subplot(2,2,4)
plot(t1,x1(:,4),'g')



%%
% X0 = x1(end,1);
% S0 = x1(end,2);
% P0 = x1(end,3);
% V = x1(end,4);
% DO = x1(end,5);
% x0 = [X0 S0 P0 V DO] ; % X S P V DO
% tspan = [5 15];
% [t2,x2] = ode15s(@(t,x) modelo_carcamo(t,x), tspan, x0); % X S P V DO
% qS_crit = qO_max./YS_ox_O ; %  [g_Glu / (g_DCW * h)]
% tasa específica de crecimiento crítica
% mu_crit = qS_crit.*YS_ox_X ; % [1/h]
% mu_set = mu_crit.*Factor_mu;
% F2 = mu_set.*Vin.*Xin.*exp(mu_set.*(t2-5))./(YS_ox_X.*Sfeed);


%F = mu_set.*Vin.*Xin.*exp(mu_set.*t)./(Yxsox.*Sfeed);


% figure(1)
% plot(t1,x1(:,1),'g',t1,x1(:,2),'r',t1,x1(:,3),'k','LineWidth',2);
% legend('Biomass','Sustrato','Ethanol')
% xlabel(' Time (h)')
% ylabel(' Concentration [g/L] ');
% ylim([0 10])
% F1 = 0;
% 
% 
% X0 = x1(end,1);
% S0 = x1(end,2);
% P0 = x1(end,3);
% V = x1(end,4);
% DO = x1(end,5);
% x0 = [X0 S0 P0 V DO] ; % X S P V DO
% tspan = [5 15];
% [t2,x2] = ode15s(@(t,x) modelo_carcamo(t,x), tspan, x0); % X S P V DO
% qS_crit = qO_max./YS_ox_O ; %  [g_Glu / (g_DCW * h)]
% tasa específica de crecimiento crítica
% mu_crit = qS_crit.*YS_ox_X ; % [1/h]
% mu_set = mu_crit.*Factor_mu;
% F2 = mu_set.*Vin.*Xin.*exp(mu_set.*(t2-5))./(YS_ox_X.*Sfeed);
% 
% figure(1)
% subplot(2,1,1)
% plot(t1,x1(:,1),'g',t1,x1(:,2),'r',t1,x1(:,3),'k','LineWidth',2);
% legend('Biomass','Sustrato','Ethanol')
% xlabel(' Time (h)')
% ylabel(' Concentration [g/L] ');
% hold on
% plot(t2,x2(:,1),'g',t2,x2(:,2),'r',t2,x2(:,3),'k','LineWidth',2);
% legend('Biomass','Sustrato','Ethanol')
% xlabel(' Time (h)')
% ylabel(' Concentration [g/L] ');
% ylim([0 50])
% 
% subplot(2,1,2);
% hold on
% yyaxis right
% plot(t2,F2,'r')
% ylabel('Flujo')
% 
% ylim([0 0.01])
% hold on
% yyaxis left
% plot(t1,x1(:,3),'k',t2,x2(:,3),'k')
% ylabel('Ethanol')
% ylim([0 2])