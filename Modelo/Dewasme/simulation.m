clear ; clc ; close all;

% Condiciones Iniciales

X0 = 0.4;   %g/L
S0 = 0.5;   %g/L
P0 = 0.8;   %g/L
O0 = 0.035; %g/L
C0 = 1.286; %g/L
V0 = 6.8;   %L

x0 = [X0 S0 P0 O0 C0 V0]';

% Parámetros de entrada al modelo

muset   = 0.29; %h-1
Sin     = 350;  %g/L

u = [muset X0 Sin V0];

% Resolver Sistema

tspan = [0 25]';
[t,x] = ode15s(@(t,x) dewasme_model(t,x,u), tspan, x0);

X = x(:,1);
S = x(:,2);
P = x(:,3);
O = x(:,4);
C = x(:,5);
V = x(:,6);

figure(1)
tiledlayout(2,1)

nexttile
plot(t,X)
hold on
grid on
plot(t,S)
plot(t,P)
legend('X','S','P','Location','west')
hold off

nexttile
plot(t,O)
hold on
grid on
plot(t,C)
legend('O','C')
hold off

figure(2)
semilogy(t,X, 'sr')

figure(3)
plot(t,V)
legend('V')
grid on
