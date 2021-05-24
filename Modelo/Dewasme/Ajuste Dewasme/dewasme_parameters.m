function dxdt = dewasme_parameters(t,x,u,k)

% Variables de Estado

dxdt    = zeros(6,1);

X       = x(1);
S       = x(2);
P       = x(3);
O       = x(4);
C       = x(5);
V       = x(6);

% Parámetros

osat    = 0.035;    % g/L
csat    = 1.286;    % g/L
kla     = 180*20;   % h-1

% Coeficiente kla tiene que ser estimado por correlaciones
% Paper: https://doi.org/10.1016/j.ijheatmasstransfer.2018.04.045
% Mientras utilicé un kla estimado de la respuesta de Carlos Araujo
% https://www.researchgate.net/post/What_is_the_common_transfer_coefficient_in_stirred_tank_bioreactor_in_both_lab_and_industrial_scale

kx1     = k(1);
kx2     = k(2);
kx3     = k(3);
ks1     = k(4);
ks2     = k(5);
kp2     = k(6);
kp3     = k(7);
ko1     = k(8);
ko2     = k(9);
ko3     = k(10);
kc1     = k(11);
kc2     = k(12);
kc3     = k(13);
muo     = k(14);
mus     = k(15);
ko      = k(16);
ks      = k(17);
kp      = k(18);
kip     = k(19);
kos     = k(20);
kop     = k(21);



% Variables manipuladas

muset   = u(1);
X0      = u(2);
Sin     = u(3);
V0      = u(4);
Fin     = muset*X0*V0*exp(muset*t)/(kx1*Sin);


% Vmax    = 20;
% if V >= 0.8*Vmax
%     Fin = 0;
% end


% Ecuaciones constitutivas

rs      = mus*(S/(S+ks));
ro      = muo*(O/(O+ko))*(kip/(kip+P));
rscrit  = ro/kos;
r1      = min(rs,rscrit)/ks1;
r2      = max(0,(rs-rscrit))/ks2;
r3      = max(0,(kos*(rscrit-rs)/kop)*(P/(P+kp)))/kp3;
D       = Fin/V;
OTR     = kla*(osat-O);
CTR     = kla*(C-csat);

% Ecuaciones Diferenciales

dxdt(1) = (kx1*r1+kx2*r2+kx3*r3)*X-D*X;     % Biomasa     dXdt
dxdt(2) = -(ks1*r1+ks2*r2)*X+D*Sin-D*S;     % Sustrato    dSdt
dxdt(3) = (kp2*r2-kp3*r3)*X-D*P;            % Producto    dPdt
dxdt(4) = -(ko1*r1+ko2*r2+ko3*r3)*X-D*O+OTR;% Oxígeno     dOdt
dxdt(5) = (kc1*r1+kc2*r2+kc3*r3)*X-D*C-CTR; % Carb. Diox. dCdt
dxdt(6) = Fin;                              % Volumen     dVdt



end