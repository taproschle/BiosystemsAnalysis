function [sys, x0,str,ts]=xu_simulink(t,x,u,flag)
% Este conjunto de lineas de codigo NO se modifica.
switch flag
   case 0
   [sys, x0,str,ts]=mdlInitializeSizes;
   case 1
   sys=mdlDerivatives(t,x,u);
   case 3
   sys=mdlOutputs(t,x,u);
   case {2,4,9}
   sys=[];
   otherwise
     error(['Unhandled flag = ',num2str(flag)]);
end
end

function [sys,x0,str,ts]=mdlInitializeSizes

% Aqui se declara el numero de ecuaciones diferenciales a integrar, de va-
% riables de entrada que ingresan al macro y de salidas que este tendra.

sizes=simsizes;
sizes.NumContStates  = 5;% Numero de ecuaciones diferenciales a integrar
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 5;% Numero de variables de salida que tendra el
%                           macro.
sizes.NumInputs      = 4; % Numero de variables de entrada que el macro
%                           aceptara.
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;

% A continuacion se presentan los valores iniciales que tendran los estados
% (las ecuaciones diferenciales a ser resueltas). Se buscara su valor en
% una primera etapa, para luego dejar fijas.

sys=simsizes(sizes); % Expresion que permite que el macro funcione bien. No
% se debe borrar.

Xin = 0.4; Sin = 0.5; Ain = 0; Oin = 2e-3; Vin = 6.8;
x0 = [Xin Sin Ain Oin Vin];

% Estas 2 lineas de codigo NO se tocan.
str=[];
ts=[0 0];
end

function sys=mdlDerivatives(~,y,u)

% Entradas y Variables de Estado
F   = u(1);         %[L/h]
N   = u(2);                         % Agitation rate [rpm]
G   = u(3);                          % Aire flow [LPM]
yo2 = u(4);                        % Fraccion gaseosa de O2

X       = y(1);         %[L]     Volumen fermentador
S       = y(2);         %[g/L]   Biomasa
A       = y(3);         %[g/L]   Glicerol
O       = y(4);         %[g/L]   1,3-Propanodiol
V       = y(5);         %[g/L]   Acido acetico

Yax = 0.667;
YS_ox_X = 0.51;
YS_of_X = 0.15;
Yoa = 1.067;
Ysa = 0.667; 
Yso = 1.067;
C_X = 0.04;
C_A = 1/30;
C_S = 1/30;
K_A = 0.05;
K_i_A = 5;
K_S = 0.05;
qm = 0.04;
qO_max = 13.4*32/1000; % g/ g h
qAc_max = 0.2;
qS_max = 1.25;

% Parámetros no ajustables:
%mu_set = 0.27;

%Xin = 0.4; Sin = 0.5; Ain = 0; Oin = 2e-3; Vin = 6.8;
Sfeed = 550;
%O_sat = 0.035; %850/1000;
K_O =  0.0001; % g o2 L-1
P = 1;

% Oxygen:

Henry =26.409;                   %  [(L*atm)/gO2] 
% Control law
alfa = 0.034;                              % [1/h], R2
beta = 1.33;                                % [1/h], R2
gamma = 0.603;                       % [1/h], R2
%delta = 0.89;                              % [-],   R3
klao2 = 5*alfa*(N^beta)*(G^gamma); % Mass transfer coefficient for O2,[1/h]
%klao2 = 180*100;


% Constitutive equations
qS = (qS_max*S/(K_S + S))*1/(1+(A/K_i_A));
qOs = min(qO_max*(O/(K_O+O))*(1/(1+(A/K_i_A))),qO_max);
qSox = min(((qOs/Yso)-(qm*YS_ox_X*C_X/C_S))/(1-YS_ox_X*C_X/C_S),qS);
qSof = max(qS - qSox,0);
qAp = (qSof - qSof*YS_of_X*C_X/C_S)*Ysa;
qAc = min(qAc_max*A/(A+K_A) , (qO_max - qOs)*Yoa/(1-Yax*C_X/C_A));
mu = (qSox - qm)*YS_ox_X + qSof*YS_of_X + qAc*Yax;
qO = qOs + (qAc - qAc*Yax*C_X/C_A)*Yoa;

%F = mu_set*Vin*Xin*exp(mu_set*t)/(YS_ox_X*Sfeed);

% Ecuaciones diferenciales.

dXdt = (mu)*X - (F/V)*X; % biomass
dSdt =  (F/V)*(Sfeed-S)- qS*X ; %sustrato
dAdt = (qAp-qAc)*X - (F/V)*A; % acetato
dOdt = klao2*(P*yo2/Henry - O)-qO*X - (F/V)*O;
dVdt = F; % volumen

sys = [dXdt dSdt dAdt dOdt dVdt];
%sys = real(sys);
end



function sys=mdlOutputs(~,y,~)

% Aqui se presenta el vector de salidas que tendra este macro, las cuales
% se ordenan para facilitar la colocacion de los bloques:

X      = y(1);     %[L] Volumen
S      = y(2);     %[g/L] Biomasa
A      = y(3);     %[g/L] Glicerol
O      = y(4);     %[g/L] 1-3 Propanodiol
V      = y(5);     %[g/L] Acido acetico


sys = [X S A O V];
sys(sys < 0) = 0;

end