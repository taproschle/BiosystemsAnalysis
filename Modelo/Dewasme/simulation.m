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
% [0.25 0.27 0.29 0.31] -> Probar
muset   = 0.29; %h-1
Sin     = 350;  %g/L

u1 = [0.25 X0 Sin V0]';
u2 = [0.27 X0 Sin V0]';
u3 = [0.29 X0 Sin V0]';
u4 = [0.31 X0 Sin V0]';
uT = [u1  u2  u3  u4];
options = odeset('NonNegative',[1, 2, 3, 4, 5, 6]);


% Resolver Sistema
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
legend('\mu = 0.25','\mu = 0.27','\mu = 0.29','\mu = 0.31','linewidth',1,'Location','northwest')
xlabel('Time [hr]')
ylabel('Biomass Concentration [g/L]')
title('Effect of \mu_{set} in Biomass and Ethanol (Dewasme Model)')
hold off

subplot(2,1,2)
for i = 1:4
    ui = uT(:,i)';
    tspan = [0 25]';
    [t,x] = ode15s(@(t,x) dewasme_model(t,x,ui), tspan, x0, options);
    P = x(:,3);
    plot(t,P,'Color',greens(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
grid on
%axis([0 25 0 1.1])
legend('\mu = 0.25','\mu = 0.27','\mu = 0.29','\mu = 0.31','linewidth',1,'Location','northwest')
xlabel('Time [hr]')
ylabel('Ethanol Concentration [g/L]')
hold off

%%
figure(1);
tspan = [0 25];
for i = 1:4
    ui = uT(:,i)';
    [t1,Y] = ode15s(@(t,x) dewasme_model(t,x,ui), tspan, x0, options);
    plot(t1, Y(:,1) ,'Color',reds(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
legend('\mu = 0.25','\mu = 0.27','\mu = 0.29','\mu = 0.31','linewidth',1,'Location','northwest')
title('Effect of \mu_{set} on Biomass (Dewasme Model)')
xlabel('Time [hr]')
ylabel('Biomass [g/l]')
grid on
figure(2);
subplot(2,1,1)
for i = 1:4
    ui = uT(:,i)';
    [t1, Y] = ode15s(@(t,x) dewasme_model(t,x,ui), tspan, x0, options);
    plot(t1, Y(:,2) ,'Color',blues(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
legend('\mu = 0.25','\mu = 0.27','\mu = 0.29','\mu = 0.31','linewidth',1,'Location','northwest')
title('Effect of \mu_{set} on Substrate (Dewasme Model)')
xlabel('Time [hr]')
ylabel('Substrate [g/l]')
grid on

subplot(2,1,2)
for i = 1:4
    ui = uT(:,i)';
    [t1, Y] = ode15s(@(t,x) dewasme_model(t,x,ui), tspan, x0, options);
    plot(t1, Y(:,3) ,'Color',greens(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
legend('\mu = 0.25','\mu = 0.27','\mu = 0.29','\mu = 0.31','linewidth',1,'Location','northwest')
title('Effect of \mu_{set} on Ethanol (Dewasme Model)')
xlabel('Time [hr]')
ylabel('Ethanol [g/l]')
grid on


%%
tspan = [0 25]';
u = [0.29 X0 Sin V0];
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

% subplot(2,1,2)
% yyaxis left
% plot(t,O,'-c','linewidth',1)
% axis([0 25 0.025 0.038])
% ylabel('O [g/L]')
% xlabel('t [hr]')
% hold on
% grid on
% yyaxis right
% ylabel('C [g/L]')
% plot(t,C,'-m','linewidth',1)
% axis([0 25 1.285 1.3])
% legend('O','C','Location','northwest')
% hold off

% figure(2)
% semilogy(t,X, 'sr')
% 
% figure(3)
% plot(t,V)
% legend('V')
% grid on
