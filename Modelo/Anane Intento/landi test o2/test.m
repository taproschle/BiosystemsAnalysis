clc,clear 
close all
% Initial conditions
X0 = 3.63; 
G0 = 0; 
E0 = 7.95; 
V0 = 2; % litros 

initials = [X0 G0 E0 V0];
% 
options = odeset("NonNegative", 1:4);
[t, y] = ode45(@landi02, [0, 25], initials', options);


%% F
Gr = 450;
t0 = 0;
Y_XG2 = 0.14./0.92;
mu_G = 0.16.*exp(-0.024.*(t-t0));
Y_XG = Y_XG2;
alpha = 0.1;
mu_set = alpha;
F0 = mu_set.*X0.*V0./(Y_XG.*Gr);
F = F0.*exp(alpha.*t);

%
X = y(:,1);
G = y(:,2);
E = y(:,3);
figure(1)
subplot(2,2,1)
plot(t, X, "LineWidth", 1)
legend('biomass')
subplot(2,2,2)
plot(t, G, "LineWidth", 1)
legend('glucose')
subplot(2,2,3)
plot(t, E, "LineWidth", 1)
legend('ethanol')
subplot(2,2,4)
plot(t, F, 'LineWidth',1)
legend('F')

% legend("Biomass", "Glucose", "Ethanol", "LineWidth", 0.7, "Location", "NorthWest")
% hold off