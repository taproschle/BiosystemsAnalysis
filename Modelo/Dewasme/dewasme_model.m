function dxdt = dewasme_model(~,x,u)

% Variables manipuladas

Fin     = u(1);     % L/h
Sin     = u(2);     % g/L

% Parámetros

osat    = 0.035;    % g/L
csat    = 1.286;    % g/L
kla     = 180;      % h-1

% Coeficiente kla tiene que ser estimado por correlaciones
% Paper: https://doi.org/10.1016/j.ijheatmasstransfer.2018.04.045
% Mientras utilicé un kla estimado de la respuesta de Carlos Araujo
% https://www.researchgate.net/post/What_is_the_common_transfer_coefficient_in_stirred_tank_bioreactor_in_both_lab_and_industrial_scale

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

% Ecuaciones constitutivas

rs      = mus*(x(2)/(x(2)+ks));
ro      = muo*(x(4)/(x(4)+ko))*(kip/(kip+x(3)));
rscrit  = ro/kos;
r1      = min(rs,rscrit)/ks1;
r2      = max(0,(rs-rscrit))/ks2;
r3      = max(0,(kos*(rscrit-rs)/kop)*(x(3)/(x(3)+kp)))/kp3;
D       = Fin/x(6);
OTR     = kla*(osat-x(4));
CTR     = kla*(x(5)-csat);

% Ecuaciones Diferenciales

dxdt    = zeros(6,1);
dxdt(1) = (kx1*r1+kx2*r2+kx3*r3)*x(1)-D*x(1);     % Biomasa     dXdt
dxdt(2) = -(ks1*r1+ks2*r2)*x(1)+D*Sin-D*x(2);     % Sustrato    dSdt
dxdt(3) = (kp2*r2-kp3*r3)*x(1)-D*x(3);            % Producto    dPdt
dxdt(4) = -(ko1*r1+ko2*r2+ko3*r3)*x(1)-D*x(4)+OTR;% Oxígeno     dOdt
dxdt(5) = (kc1*r1+kc2*r2+kc3*r3)*x(1)-D*x(5)-CTR; % Carb. Diox. dCdt
dxdt(6) = Fin;                                    % Volumen     dVdt

end