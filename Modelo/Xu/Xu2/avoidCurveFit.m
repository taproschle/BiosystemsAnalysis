function out = avoidCurveFit(k, t)

times = t;
y0 = [5 0.04 0 0.004 0.3];
% Parámetros no ajustables:
mu_set = 0.13;
klao2 = 180*100;
Xin = 5; Sin = 0.04; Ain = 0; Oin = 4e-3; Vin = 0.3;
Sfeed = 550;
O_sat = 0.035; %850/1000;
K_O =  0.0001; % g o2 L-1

v = [mu_set klao2 Vin Xin Sfeed O_sat K_O];


[~, y] = ode45(@(t,x) xu_model(t,x,v,k), times, y0);

X = y(:,1);
S = y(:,2);

out = [X S];
end