
% Orden parametros [X S A DOT DOTa V]
initial = [0.17 4.94 0.0129 99 99 2];
time_span = [0 35];

[t,x] = ode45(@anane_model, time_span, initial);

% Ploteo de resultados

X = x(:,1);
S = x(:,2);
A = x(:,3);
DOT = x(:,4);
DOTa = x(:,5);
V = x(:,6);

%% Subplots

%subplot(1,3,1)
plot(t, X, "-", "LineWidth", 2)
hold on
plot(t, S, "-", "LineWidth", 2)
plot(t, A, "-", "LineWidth", 2)
legend("Biomass", "Sustrate", "Acetate", "Location", "northwest", "LineWidth",0.7)
xlabel("Time [h]")
ylabel("Concentration [g/L]")
title("Evolution of compounds through a fed-batch model")
hold off

%% A
subplot(1,3,2)
plot(t, DOT)
hold on
plot(t, DOTa)
hold off
subplot(1,3,3)
plot(t, V)

