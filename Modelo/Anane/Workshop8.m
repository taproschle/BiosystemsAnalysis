%% Variables de holguras 

z = @(x) 200*exp(-0.02.*sqrt(x(1).^2+x(2).^2)) + ...
    5*exp(cos(3*x(1))+sin(3*x(2)));

z2 = @(x,y) 200*exp(-0.02.*sqrt(x.^2+y.^2)) + ...
    5*exp(cos(3*x)+sin(3*y));
% For case a)

x1 = -32:32;
x2 = -32:32;

[X,Y] = meshgrid(x1,x2);


surf(X, Y, z2(X,Y))
hold on
plot3(-32, 30.97, 83.65, "dr", "MarkerSize", 10, "MarkerFaceColor", "r")
plot3(32, 30.97, 83.65, "dr", "MarkerSize", 10, "MarkerFaceColor", "r")
plot3(3.17, -32, 105.81, "dr", "MarkerSize", 10, "MarkerFaceColor", "r")
plot3(-32, -29.92, 84.9, "dr", "MarkerSize", 10, "MarkerFaceColor", "r")
plot3(32, -29.92, 84.9, "dr", "MarkerSize", 10, "MarkerFaceColor", "r")
xlabel("Eje x")
ylabel("Eje y")
zlabel("Valor función")
hold off
%% For case b)

x1 = 15:30;
x2 = 15:30;

x0 = [0, 0];
%x0 = [-30, 30];
%x0 = [30, 30];
%x0 = [-30, -30];
%x0 = [30, -30];
[x_sal, fval, flag] = fmincon(z, x0, [], [], [], [], [-32 -32], [32 32]);

fprintf("Vector salida %f \n", x_sal)
fprintf("funcion objetivo %f \n", fval)


%% Caso maximización
x1 = -32:0.4:32;
x2 = -32:0.4:32;

[X,Y] = meshgrid(x1,x2);

x0 = [0, 0];
%x0 = [-30, 30];
%x0 = [30, 30];
%x0 = [-30, -30];
%x0 = [30, -30];
z_max = @(x) -200*exp(-0.02.*sqrt(x(1).^2+x(2).^2)) - ...
    5*exp(cos(3*x(1))+sin(3*x(2)));
[x_sal, fval, flag] = fmincon(z_max, x0, [], [], [], [], [-32 -32], [32 32]);
disp(fval)
disp(x_sal)
% Solo ploteo


surf(X, Y, z2(X,Y));

hold on
plot3(0, 0.5117, 234.88, "dr", "MarkerSize", 10, "MarkerFaceColor", "r")
plot3(-6.2, 21.5, 164.81, "dr", "MarkerSize", 10, "MarkerFaceColor", "r")
plot3(6.2, 21.46, 164.81, "dr", "MarkerSize", 10, "MarkerFaceColor", "r")
plot3(20.93, -22.51, 145.08, "dr", "MarkerSize", 10, "MarkerFaceColor", "r")
plot3(-20.93, -22.51, 145.08, "dr", "MarkerSize", 10, "MarkerFaceColor", "r")
xlabel("Eje x")
ylabel("Eje y")
zlabel("Valor función")
hold off




