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
sizes.NumOutputs     = 8;% Numero de variables de salida que tendra el
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

load vXu.mat v

X0 = v(2);
V0 = v(3);
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;

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

X       = y(1);         %[L]     Volumen fermentador
S       = y(2);         %[g/L]   Biomasa
E       = y(3);         %[g/L]   Glicerol
O       = y(4);         %[g/L]   1,3-Propanodiol
V       = y(5);         %[g/L]   Acido acetico

load vXu.mat v
load kXu.mat k

% Parámetros no ajustables:
Sin       = v(4);
Ko        = v(7);
Kio       = v(8);

Ks      = k(1);
qSmax   = k(2);
Ysoxx   = k(3);
qm      = k(4);
Yso     = k(5);
Kie     = k(6);
Yse     = k(7);
Kec     = k(8);
qEcmax   = k(9);
Ysofx   = k(10);
Yeo     = k(11);
Yex     = k(12);
qOmax   = k(13);

P = 1;
Henry =26.409;                   %  [(L*atm)/gO2] 
% Control law
alfa = 0.034;                              % [1/h], R2
beta = 1.33;                                % [1/h], R2
gamma = 0.603;                       % [1/h], R2
%delta = 0.89;                              % [-],   R3
klao2 = 5*alfa*(N^beta)*(G^gamma); % Mass transfer coefficient for O2,[1/h]

% Constitutive equations
qS      = (qSmax*S/(Ks+S))*1/(1+(E/Kie));
qOs     = qOmax*(O/(Ko+O))*(1/(1+(E/Kio)));
qSox    = min(((qOs/Yso)-(qm*Ysoxx))/(1-Ysoxx),qS);
qSof    = max(qS-qSox,0);
qEp     = (qSof-qSof*Ysofx)*Yse;
qEc     = min(qEcmax*E/(E+Kec),(qOmax-qOs)*Yeo/(1-Yex));
mu      = (qSox-qm)*Ysoxx+qSof*Ysofx+qEc*Yex;
qO      = qOs+(qEc-qEc*Yex)*Yeo;
D       = F/V;
% Ecuaciones diferenciales.

dXdt = X*(mu-D);                 % dXdt
dSdt = D*(Sin-S)-qS*X;           % dSdt
dEdt = (qEp-qEc)*X-D*E;          % dEdt
dOdt = klao2*(P*yo2/Henry-O)-qO*X-D*O;  % dOdt
dVdt = F; % Volume                 % dVdt

sys = [dXdt dSdt dEdt dOdt dVdt];
end

function sys=mdlOutputs(~,y,~)
% Aqui se presenta el vector de salidas que tendra este macro, las cuales
% se ordenan para facilitar la colocacion de los bloques:

X      = y(1);     %[L] Volumen
S      = y(2);     %[g/L] Biomasa
E      = y(3);     %[g/L] Glicerol
O      = y(4);     %[g/L] 1-3 Propanodiol
V      = y(5);     %[g/L] Acido acetico

% Parámetros no ajustables:

load vXu.mat v
load kXu.mat k

Ko        = v(7);
Kio       = v(8);

% Parámetros ajustables
Ks      = k(1);
qSmax   = k(2);
Ysoxx   = k(3);
qm      = k(4);
Yso     = k(5);
Kie     = k(6);
Kec     = k(8);
qEcmax  = k(9);
Ysofx   = k(10);
Yeo     = k(11);
Yex     = k(12);
qOmax   = k(13);

% Constitutive equations
qS      = (qSmax*S/(Ks+S))*1/(1+(E/Kie));
qOs     = qOmax*(O/(Ko+O))*(1/(1+(E/Kio)));
qSox    = min(((qOs/Yso)-(qm*Ysoxx))/(1-Ysoxx),qS);
qSof    = max(qS-qSox,0);
qEc     = min(qEcmax*E/(E+Kec),(qOmax-qOs)*Yeo/(1-Yex));
mu      = (qSox-qm)*Ysoxx+qSof*Ysofx+qEc*Yex;
qO      = qOs+(qEc-qEc*Yex)*Yeo;
qScrit  = ((qOs/Yso)-(qm*Ysoxx))/(1-Ysoxx);
mu_crit = qScrit*Ysoxx;
Scrit   = Ks*qO/(Yso*qSmax - qO);

sys = [X S E O V mu mu_crit Scrit];
sys(sys < 0) = 0;
end