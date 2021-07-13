function [sys, x0,str,ts]=an_simulink(t,x,u,flag)
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
sizes.NumInputs      = 19; % Numero de variables de entrada que el macro
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

X       = y(1);         %[L]     Volumen fermentador
S       = y(2);         %[g/L]   Biomasa
E       = y(3);         %[g/L]   Glicerol
O       = y(4);         %[g/L]   1,3-Propanodiol
V       = y(5);         %[g/L]   Acido acetico

load vAn.mat v
% load kAn.mat k

% Parámetros no ajustables:
Sin     = v(4);
Ko      = v(7);


Ks      = u(5);
qm      = u(6);
Yso     = u(7);
Kie     = u(8);
qEpmax  = u(9);
Kep     = u(10);
Yse     = u(11);
Kis     = u(12);
Yeo     = u(13);
qSmax   = u(14);
qEcmax  = u(15);
Kec     = u(16);
Ysoxx   = u(17);
Ysofx   = u(18);
Yex     = u(19);

P = 1;
Henry =26.409;                   %  [(L*atm)/gO2] 
% Control law
alfa = 0.034;                              % [1/h], R2
beta = 1.33;                                % [1/h], R2
gamma = 0.603;                       % [1/h], R2
%delta = 0.89;                              % [-],   R3
klao2 = 5*alfa*(N^beta)*(G^gamma); % Mass transfer coefficient for O2,[1/h]

% Constitutive equations
qS      = (qSmax/(1+(E/Kie)))*(S/(S+Ks));
qSof    = qEpmax*(qS/(qS+Kep));
qEp      = qSof*Yse;
qSox    = (qS-qSof)*(O/(O+Ko));
qEc     = (qEcmax/(1+(qS/Kis)))*(E/(E+Kec));
qE      = qEp-qEc;
mu      = (qSox-qm)*Ysoxx+qSof*Ysofx+qEc*Yex;
qO      = Yso*(qSox-qm)+qEc*Yeo;
D       = F/V;

% Ecuaciones diferenciales.

dXdt = X*(mu-D);                 % dXdt
dSdt = D*(Sin-S)-qS*X;           % dSdt
dEdt = qE*X-D*E;                 % dEdt
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

% load vAn.mat v
% load kAn.mat k
% 
% % Parámetros no ajustables:
% Ko      = v(7);
% 
% % Fixed parameters (all cases)
% Ks      = 0.0450;
% qm      = 0.0067;
% Yso     = 0.4684;
% Kie     = 1.6129;
% qEpmax  = 0.6350;
% Kep     = 0.5034;
% Kis     = 2.2729;
% Yeo     = 0.8163;
% 
% % Adjusted parameters (all cases)
% qSmax   = k(1);
% qEcmax  = k(2);
% Kec     = k(3);
% Ysoxx   = k(4);
% Ysofx   = k(5);
% Yex     = k(6);
% 
% % Constitutive equations
% qS      = (qSmax/(1+(E/Kie)))*(S/(S+Ks));
% qSof    = qEpmax*(qS/(qS+Kep));
% qSox    = (qS-qSof)*(O/(O+Ko));
% qEc     = (qEcmax/(1+(qS/Kis)))*(E/(E+Kec));
% mu      = (qSox-qm)*Ysoxx+qSof*Ysofx+qEc*Yex;
% qO      = Yso*(qSox-qm)+qEc*Yeo;
% qscrit  = qS*O/(O+Ko);
% mucrit  = qscrit*Ysoxx;
% scrit   = Ks*qO/(Yso*qSmax-qO);

sys = [X S E O V];
sys(sys < 0) = 0;
end