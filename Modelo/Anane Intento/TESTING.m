clc
% Initial conditions
X0 = 3.63; 
G0 = 0; 
E0 = 7.95; 
V0 = 2; % litros 

initials = [X0 G0 E0 V0];
% 
options = odeset("NonNegative", 1:4);
[t, y] = ode45(@Landi, [0, 25], initials', options);

%
X = y(:,1);
G = y(:,2);
E = y(:,3);

plot(t, X, "LineWidth", 1)
hold on
plot(t, G, "LineWidth", 1)
plot(t, E, "LineWidth", 1)

legend("Biomass", "Glucose", "Ethanol", "LineWidth", 0.7, "Location", "NorthWest")
hold off