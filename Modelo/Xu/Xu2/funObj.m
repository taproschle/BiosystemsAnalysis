function f = funObj(k)

tabla = readtable('datos_troles.csv');
tabla = table2array(tabla);

texp = tabla(:,1);
yexp = tabla(:,2:6);


mu_set = 0.13;
klao2 = 180*100;
Xin = 5; Vin = 0.3;
Sfeed = 550;
O_sat = 0.035; %850/1000;
K_O =  0.0001; % g o2 L-1

v = [mu_set klao2 Vin Xin Sfeed O_sat K_O];
y0 = [5 0.04 0 0.004 0.3];


tspan = texp;
fun = @(t,y) xu_model(t,y,v,k);
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'NonNegative',[1,2,3,4,5]);
disp('Integrando')
disp(randi(100))
[~ , Y] = ode15s(fun,tspan,y0,options);

pesos = [1 1 1 1 1];

if ~ isequal(size(Y),size(yexp))
    n = 1e10*ones(1,5);
    f = dot(pesos,n);
else
    nX = norm(Y(:,1) - yexp(:,1));
    nS = norm(Y(:,2) - yexp(:,2));
    nA = norm(Y(:,3) - yexp(:,3));
    nO = norm(Y(:,4) - yexp(:,4));
    nV = norm(Y(:,5) - yexp(:,5));
    n = [nX nS nA nO nV].^2;
    f = dot(pesos,n);
end
end