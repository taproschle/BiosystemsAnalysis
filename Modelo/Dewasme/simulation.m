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
% [0.25 0.27 0.29 0.31] -> Probar
muset   = 0.29; %h-1
Sin     = 350;  %g/L

u1 = [0.25 X0 Sin V0]';
u2 = [0.27 X0 Sin V0]';
u3 = [0.29 X0 Sin V0]';
u4 = [0.231 X0 Sin V0]';
uT = [u1  u2  u3  u4];
options = odeset('NonNegative',[1, 2, 3, 4, 5, 6]);

% Resolver Sistema

c1 = [0.99607 0.850984 0.46274]';
c2 = [0.99215 0.552941 0.23529]';
c3 = [0.98823 0.305882 0.16470]';
c4 = [0.89019 0.101960 0.10980]';
reds = [c1 c2 c3 c4];
b1 = [0.61967 0.792157 0.88235]';
b2 = [0.41960 0.682353 0.83921]';
b3 = [0.25882 0.572549 0.77647]';
b4 = [0.12941 0.443137 0.70980]';
blues = [b1 b2 b3 b4];
g1 = [0.63137 0.850980 0.60784]';
g2 = [0.45902 0.768627 0.46274]';
g3 = [0.25490 0.670588 0.36470]';
g4 = [0.00000 0.352941 0.19607]';
greens = [g1 g2 g3 g4];

forms = ["--" , "-." , "-" , "--"];

figure(1);
subplot(2,1,1)
for i = 1:4
    ui = uT(:,i)';
    tspan = [0 25]';
    [t,x] = ode15s(@(t,x) dewasme_model(t,x,ui), tspan, x0, options);
    X = x(:,1);
    plot(t,X,'Color',reds(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
grid on
legend('\mu = 0.25','\mu = 0.27','\mu = 0.29','\mu = 0.31','linewidth',1,'Location','west')
xlabel('Time [hr]')
ylabel('Biomass Concentration [g/L]')
title('Effect of \mu_{set} in Biomass and Ethanol Concentrations')
hold off

subplot(2,1,2)
for i = 1:4
    ui = uT(:,i)';
    tspan = [0 25]';
    [t,x] = ode15s(@(t,x) dewasme_model(t,x,ui), tspan, x0, options);
    P = x(:,3);
    plot(t,x(:,3),'Color',greens(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
grid on
%axis([0 25 0 1.1])
legend('\mu = 0.25','\mu = 0.27','\mu = 0.29','\mu = 0.31','linewidth',1)
xlabel('Time [hr]')
ylabel('Ethanol Concentration [g/L]')
hold off


%%
tspan = [0 25]';
u = [0.31 X0 Sin V0];
[t,x] = ode15s(@(t,x) dewasme_model(t,x,u), tspan, x0, options);

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
subplot(2,1,1)
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
legend('X','S','P','Location','west')
title('State Variables in Dewasme Model')
xlabel('Time [hr]')
ylabel('S and P [gr/L]')
%axis([0 25 0 1.2])
hold off

subplot(2,1,2)
yyaxis left
plot(t,O,'-c','linewidth',1)
axis([0 25 0.025 0.038])
ylabel('O [g/L]')
xlabel('t [hr]')
hold on
grid on
yyaxis right
ylabel('C [g/L]')
plot(t,C,'-m','linewidth',1)
axis([0 25 1.285 1.3])
legend('O','C','Location','west')
hold off

% figure(2)
% semilogy(t,X, 'sr')
% 
% figure(3)
% plot(t,V)
% legend('V')
% grid on
