function [sys, x0,str,ts]=dew_simulink(t,x,u,flag)
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
sizes.NumOutputs     = 8;% Numero de variables de salida que tendra el
%                           macro.
sizes.NumInputs      = 16; % Numero de variables de entrada que el macro
%                           aceptara.
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;

% A continuacion se presentan los valores iniciales que tendran los estados
% (las ecuaciones diferenciales a ser resueltas). Se buscara su valor en
% una primera etapa, para luego dejar fijas.

sys=simsizes(sizes); % Expresion que permite que el macro funcione bien. No
% se debe borrar.

X0 = 4.125;
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;
V0 = 0.3;
x0 = [X0 S0 E0 O0 V0];

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
% Parámetros 
Ks      = u(5);
qSmax   = u(6);
Ysoxx   = u(7);
Yso     = u(8);
Kio     = u(9);
Yse     = u(10);
Kec     = u(11);
Ysofx   = u(12);
Yeo     = u(13);
Yex     = u(14);
qOmax   = u(15);
Yosof   = u(16);

X       = y(1);         %[L]     Volumen fermentador
S       = y(2);         %[g/L]   Biomasa
E       = y(3);         %[g/L]   Glicerol
O       = y(4);         %[g/L]   1,3-Propanodiol
V       = y(5);         %[g/L]   Acido acetico

load vDew.mat v
load kDew.mat k

% Parámetros no ajustables:
Sin     = v(4);
Ko      = v(7);

P = 1;
Henry =26.409;                   %  [(L*atm)/gO2] 
% Control law
alfa = 0.034;                              % [1/h], R2
beta = 1.33;                                % [1/h], R2
gamma = 0.603;                       % [1/h], R2
%delta = 0.89;                              % [-],   R3
klao2 = 5*alfa*(N^beta)*(G^gamma); % Mass transfer coefficient for O2,[1/h]

% Constitutive equations
qS      = qSmax*S/(S+Ks);
qO      = (qOmax*O/(O+Ko))*(Kio/(Kio+E));
qScrit  = qO/Yso;
qSox    = min(qS,qScrit);
qSof    = max(0,qS-qScrit);
qE      = max(0,(Yso/Yeo)*(qScrit-qS)*(E/(E+Kec)));
mu      = Ysoxx*qSox + Ysofx*qSof + Yex*qE;
D       = F/V;

% Ecuaciones diferenciales.

dXdt = X*(mu-D);                 % dXdt
dSdt = D*(Sin-S)-(qSox-qSof)*X;  % dSdt
dEdt = (Yse*qSof-qE)*X-D*E;      % dEdt
dOdt = klao2*(P*yo2/Henry-O)-D*O-...    % dOdt
    (Yso*qSox+Yosof*qSof+Yeo*qE)*X;
dVdt = F; % Volume                 % dVdt

sys = [dXdt dSdt dEdt dOdt dVdt];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sys=mdlOutputs(~,y,~)
% Aqui se presenta el vector de salidas que tendra este macro, las cuales
% se ordenan para facilitar la colocacion de los bloques:

X      = y(1);     %[L] Volumen
S      = y(2);     %[g/L] Biomasa
E      = y(3);     %[g/L] Glicerol
O      = y(4);     %[g/L] 1-3 Propanodiol
V      = y(5);     %[g/L] Acido acetico

sys = [X S E O V];
sys(sys < 0) = 0;

% Constitutive equations

load vDew.mat v
load kDew.mat k

% Parámetros no ajustables:
Ko      = v(7);
Ks      = k(1);
qSmax   = k(2);
Ysoxx   = k(3);
Yso     = k(4);
Kio     = k(5);
Kec     = k(7);
Ysofx   = k(8);
Yeo     = k(9);
Yex     = k(10);
qOmax   = 0.3; %k(11);

% Constitutive equations:
qS      = qSmax*S/(S+Ks);
qO      = (qOmax*O/(O+Ko))*(Kio/(Kio+E));
qScrit  = qO/Yso;
qSox    = min(qS,qScrit);
qSof    = max(0,qS-qScrit);
qE      = max(0,(Yso/Yeo)*(qScrit-qS)*(E/(E+Kec)));
mu      = Ysoxx*qSox + Ysofx*qSof + Yex*qE;
mu_crit = qScrit*Ysoxx;
Scrit   = Ks*qO/(Yso*qSmax - qO);

sys = [sys mu mu_crit Scrit];
sys(sys < 0) = 0;
end