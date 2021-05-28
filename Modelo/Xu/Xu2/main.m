% Main Xu 
clc ; clear
close all

% Parámetros ajustables:
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

% Parámetros no ajustables:
mu_set = 0.27;
klao2 = 180*100;
Xin = 0.4; Sin = 0.5; Ain = 0; Oin = 2e-3; Vin = 6.8;
Sfeed = 550;
O_sat = 0.035; %850/1000;
K_O =  0.0001; % g o2 L-1

% Lista de parámetros ajustables entregados al modelo
k = [YS_ox_X YS_of_X Yax Yoa Ysa Yso C_X C_A C_S qm ...
    qS_max qAc_max qO_max K_S K_i_A K_A];

% Lista de parámetros no ajustables:
v = [mu_set klao2 Vin Xin Sfeed O_sat K_O];

% Resolución ODE:

tspan = [0 25];
x0 = [Xin Sin Ain Oin Vin];
fun = @(t,y) xu_model(t,y,v,k);
[T,C] = ode23s(fun, tspan, x0);

% Gráfico:
figure(1)
subplot(2,2,1)
yyaxis left
plot(T,C(:,2),'g',"LineWidth",2);
ylabel('Glucose [g/L]')
yyaxis right
plot(T,C(:,1),'r',"LineWidth",2);
ylabel('Biomass [g/L]');
legend('Glucose','Biomass');
xlabel('Time [h]')
xlim([0 25])

subplot(3,2,2)
plot(T,C(:,3),'k','LineWidth',2);
legend('Acetate')
xlabel('Time [h]')
ylabel('Acetate [g/L]');
xlim([0 25])

subplot(2,2,3)
plot(T,C(:,5),'g','LineWidth',2);
legend('Volumen [L]')
xlabel('Time [h]')
ylabel('Volume [L]');
xlim([0 25])


subplot(2,2,4)
plot(T,C(:,4),'g','LineWidth',2)
legend('Oxygen [g/L]')
xlabel('Time [h]')
ylabel('Oxygen [g/L]')
xlim([0 25])


