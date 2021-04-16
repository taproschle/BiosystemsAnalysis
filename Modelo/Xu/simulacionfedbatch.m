clc,clear
close all
global Yas Yoa Yos Yxa Yxsof Yxsox Ka KiS Ks Ko qAcmax qOmax qSmax Osat kla Vin Xin mumax



% parametros
Yas = 0.23; %  g A /g S stochiometric cte CARCAMO
%Yas = 0.667; % XU
Yoa = 20.66/1000; % g o2 / g A CARCAMO
%Yoa = 1.067; % XU 

Yos = 401/1000; % g o2 / g S CARCAMO
%Yos = 1.067; % XU
Yxa = 14.4; % g DCW / gIP CARCAMO CREO Q ESTA MALO
%Yxa = 0.4; % XU
Yxsof = 0.7; %CARCAMO
%Yxsof = 0.15; % XU
Yxsox = 0.3; % CARCAMO
%Yxsox = 0.51; %XU


Ca = 1/30; % mol Carbon / g acetato
Cs = 1/30; % mol Carbon /g azucar
Cx = 0.04; % mol Carbon / g biomasa


%KiO = 4; % g/L fitting
Ka = 0.1; % CARCAMO
%Ka = 0.05; %XU
KiS = 0.45; % g A /L fitting CARCAMO
%KiS = 5; % XU
Ks = 8; % g S/L CARCAMO
%Ks = 0.05; % XU
Ko = 0.0045; %g O/L CARCAMO
%Ko = 0.0001; % DEWASME g/L
qAcmax = 6; % g/g hr CARCAMO
%qAcmax = 0.06; % XU
qm = 0.04; % g/ g hr
qOmax = 83.5/1000; % g/g hr CARCAMO
%qOmax = 15.6*32/1000; % XU
qSmax = 4.9; % g/g hr CARCAMO
%qSmax = 1.3; % Xu
muc = 0.14;
Osat = 0.035;
kla = 202*20; % [h-1]
tspan = [0 15];
mumax = 0.55;
Vin = 6.8;
Xin = 0.1;

mu_set = 0.2;
Sfeed = 520;
u =[mu_set Xin Vin Sfeed]; % mu_Set XIN Vin Sin


x0 = [7.6 0 Xin Vin 0.035]; % S Et X V O

[t,x] = ode15s(@(t,x) fedbatch(t,x,u), tspan, x0,u);

F = mu_set.*Vin.*Xin.*exp(mu_set.*t)./(Yxsox.*Sfeed);


figure(1)
plot(t,x(:,1),'g',t,x(:,3),'r',t,x(:,2),'k','LineWidth',2);
legend('Glucose','Biomass','Ethanol')
xlabel(' Time (h)')
ylabel(' Concentration [g/L] ');
ylim([0 30])
%%
% subplot(2,2,2)
% plot(t,x(:,4),'b','LineWidth',2)
% legend('V [L]');
% xlabel('Time (h)')
% ylabel(' Volume [L]')

% %subplot(2,2,3)
% 
% %yyaxis left
% %plot(t,x(:,2),'b',"LineWidth",2);
% %ylim([0 2]);
% %ylabel('Ethanol [g/L]')
% %yyaxis right
% %plot(t,F(:,1),'r',"LineWidth",2);
% %ylabel('F [L/h]');
% %legend('Ethanol','F(t)');
% %!ylim([0 0.1]);
% %xlabel('Time [h]')
% 
% %%%


%%%
% 
% subplot(2,2,4)
% plot(t,F(:,1),'k',"LineWidth",2);
% legend('F(t) [L/h]')

%subplot(3,3,5)
%plot(t,mu,'g',"LineWidth",2)
%legend('mu [h-1]')
