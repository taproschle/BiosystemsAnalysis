clear ; clc ; close all;
colors
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

u1 = [muset X0 Sin V0]';

% Parámetros del modelo:
kx1     = 0.49;
kx2     = 0.05;
kx3     = 0.72;
ks1     = 1;
ks2     = 1;
kp2     = 0.48;
kp3     = 1;
ko1     = 0.3968;
ko2     = 0;
ko3     = 1.104;
kc1     = 0.5897;
kc2     = 0.4621;
kc3     = 0.6249;
muo     = 0.256;
mus     = 3.5;
ko      = 0.0001;
ks      = 0.1;
kp      = 0.1;
kip     = 10;
kos     = ko1;
kop     = ko3;

k = [kx1 kx2 kx3 ks1 ks2 kp2 kp3 ko1 ko2 ko3 kc1 kc2 kc3 muo mus ko ks ...
    kp kip kos kop];

options = odeset('NonNegative',[1, 2, 3, 4, 5, 6]);

tspan = [0 25]';
u = [0.29 X0 Sin V0];

[t,x] = ode15s(@(t,x) dewasme_parameters(t,x,u,k), tspan, x0, options);

X = x(:,1);
S = x(:,2);
P = x(:,3);
O = x(:,4);
C = x(:,5);
V = x(:,6);

fig = figure(2);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
plot(t,X,'-r','linewidth',1)
hold on
ylabel('Biomass [gr/L]')
xlabel('Time [hr]')


grid on
yyaxis right
plot(t,S,'-b','linewidth',1)
hold on
plot(t,P,'-g','linewidth',1)
legend('X','S','P','Location','northwest')
title('State Variables in Dewasme Model')
xlabel('Time [hr]')
ylabel('S and P [gr/L]')
%axis([0 25 0 1.2])
hold off
