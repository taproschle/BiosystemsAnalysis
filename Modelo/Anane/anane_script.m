clear ; clc ; close all
% Orden parametros [X S A DOT DOTa V]
initial = [0.5 0.5 0.0129 99 99 2];
time_span = [0 25];

% Parametros de entrada
X0 = 0.4;
V0 = 6.8;
Sin = 350;
muset = 0.22;

u = [X0 V0 Sin muset];

options = odeset('NonNegative',[1 2 3 4 5 6]);

[t,x] = ode15s(@(t,x) anane_model(t,x,u), time_span, initial, options);

% Resultados

X = x(:,1);
S = x(:,2);
A = x(:,3);
DOT = x(:,4);
DOTa = x(:,5);
V = x(:,6); 

% Subplots

%subplot(1,3,1)
plot(t, X, "-", "LineWidth", 2)
grid on
hold on
plot(t, S, "-", "LineWidth", 2)
plot(t, A, "-", "LineWidth", 2)
legend("Biomass", "Sustrate", "Acetate", "Location", "northwest", "LineWidth",0.7)
xlabel("Time [h]")
ylabel("Concentration [g/L]")
title("Evolution of compounds through a fed-batch model")
% ylim([-1, 50])
%xline(15.7, "--", "Color", [0.5, 0.5, 0.5],'HandleVisibility','off')
xline(7, "--", "Color", [0.5, 0.5, 0.5],'HandleVisibility','off')
hold off

%% Ploteo de DOT y el volumen
subplot(1,2,1)
plot(t, DOT)
hold on
plot(t, DOTa)
hold off
subplot(1,2,2)
plot(t, V, "LineWidth", 1.2)
hold on
plot(T, Y(:,1), "LineWidth", 1.2)
legend("Mio", "anane")
hold off