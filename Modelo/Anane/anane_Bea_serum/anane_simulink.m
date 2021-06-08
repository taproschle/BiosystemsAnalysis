function [sys, x0,str,ts]=anane_simulink(t,x,u,flag)
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
sizes.NumOutputs     = 6;% Numero de variables de salida que tendra el
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

X       = y(1);         %[L]     X
S       = y(2);         %[g/L]   S
A       = y(3);         %[g/L]  A
O       = y(4);         %[g/L] O2
V       = y(5);         %[g/L]   V


%PARAMETROS
   %[Kap     Ksa     Ko      Ks     Kia      Kis    pAmax   qAmax    qm     qSmax    Yas     Yoa     Yxa    Yem     Yos    Yxsof] 
Par = [0.5088	0.0128	0.1	0.0381	1.2602	1.8383	0.2286	0.1148	0.0133	0.635	0.8938	0.5221	0.5794	0.5321	1.5722	0.229];
Kap     = Par(1);   %monod-type saturation constant, intracellular acetate production
Ksa     = Par(2);   %affinity constant, acetate consumption (0.05 g/L)
Ko      = Par(3);   %affinity constant, oxygen consumption (molO./L)
Ks      = Par(4);   %affinity constant, substrate consumption (0.05 gglu./L)
Kia     = Par(5);   %inhibition constant, inhib. of glucose uptake by acetate
Kis     = Par(6);   %inhibition constant, inhib. of ace.uptake by glucose(4g/L)
pAmax   = Par(7);   %max spec acetate production rate (0.15 g_ace/(gx.h)))
qAmax   = Par(8);   %max spec acetate consumption rate (0.15 g_ace/(gx.h)))
qm      = Par(9);   %spec maintenance coefficient (0.04 g_glu/(gx.h))
qSmax   = Par(10);  %max spec glucose uptake rate (1.3 g_glu./gx.h)
Yas     = Par(11);  %yield of acetate on substrate (g.ace/g.glu)
Yoa     = Par(12);  %yield of oxygen on acetate g/g)
Yxa     = Par(13);  %yield of biomass on acetate
Yem     = Par(14);  %yield exclusive maintenance (Yem = Yxs-qm)
Yos     = Par(15);  %yield of oxygen on glucose (g/g)
Yxsof   = Par(16);  %biomass yield from the overflow route
% qOmax = 13.4*32/1000; % g g-1 h-1


% Parámetros no ajustables:
Sfeed = 550; % Sfeed del ley de control
Ko =  0.0001;% g o2 L-1 monod O0
P = 1; %presion sistema

% Oxygen:

Henry =26.409;                   %  [(L*atm)/gO2] 
% Control law
alfa = 0.034;                              % [1/h], R2
beta = 1.33;                                % [1/h], R2
gamma = 0.603;                       % [1/h], R2
%delta = 0.89;                              % [-],   R3
klao2 = 5*alfa*(N^beta)*(G^gamma); % Mass transfer coefficient for O2,[1/h]

% Constitutive equations
qS      = (qSmax/(1+(A/Kia)))*(S/(S+Ks)); % Check
qSof    = pAmax*(qS/(qS+Kap));             % Check
pA      = qSof*Yas;                        % Check
qSox    = (qS-qSof)*(O/(O+Ko));        % Check      
%qSan    = (qSox-qm)*Yem*Cx/Cs;             % 
qsA     = (qAmax/(1+(qS/Kis)))*(A/(A+Ksa));
qA      = pA-qsA;                            
mu      = (qSox-qm)*Yem + qsA*Yxa + qSof*Yxsof;
qO      = Yos*(qSox-qm)+qsA*Yoa;
% Ecuaciones diferenciales.

dXdt = (mu)*X - (F/V)*X; % biomass
dSdt =  (F/V)*(Sfeed-S)- qS*X ; %sustrato
dAdt = qA*X - (F/V)*A; % acetato
dOdt = klao2*(P*yo2/Henry - O)-qO*X - (F/V)*O;
dVdt = F; % Volume
%dCdt = -klaco2*(CO2-CO2s) + rc*X - CO2*D;  % CO2 balance [g/L/h]

sys = [dXdt dSdt dAdt dOdt dVdt];
end

function sys=mdlOutputs(~,y,~)
% Aqui se presenta el vector de salidas que tendra este macro, las cuales
% se ordenan para facilitar la colocacion de los bloques:


X      = y(1);     %[L] Volumen
S      = y(2);     %[g/L] Biomasa
A      = y(3);     %[g/L] Glicerol
O      = y(4);     %[g/L] 1-3 Propanodiol
V      = y(5);     %[g/L] Acido acetico
Par = [0.5088	0.0128	0.0001	0.0381	1.2602	1.8383	0.2286	0.1148	0.0133	0.635	0.8938	0.5221	0.5794	0.5321	1.5722	0.229];

Kap     = Par(1);   %monod-type saturation constant, intracellular acetate production
Ksa     = Par(2);   %affinity constant, acetate consumption (0.05 g/L)
Ko      = Par(3);   %affinity constant, oxygen consumption (molO./L)
Ks      = Par(4);   %affinity constant, substrate consumption (0.05 gglu./L)
Kia     = Par(5);   %inhibition constant, inhib. of glucose uptake by acetate
Kis     = Par(6);   %inhibition constant, inhib. of ace.uptake by glucose(4g/L)
pAmax   = Par(7);   %max spec acetate production rate (0.15 g_ace/(gx.h)))
qAmax   = Par(8);   %max spec acetate consumption rate (0.15 g_ace/(gx.h)))
qm      = Par(9);   %spec maintenance coefficient (0.04 g_glu/(gx.h))
qSmax   = Par(10);  %max spec glucose uptake rate (1.3 g_glu./gx.h)
Yas     = Par(11);  %yield of acetate on substrate (g.ace/g.glu)
Yoa     = Par(12);  %yield of oxygen on acetate g/g)
Yxa     = Par(13);  %yield of biomass on acetate
Yem     = Par(14);  %yield exclusive maintenance (Yem = Yxs-qm)
Yos     = Par(15);  %yield of oxygen on glucose (g/g)
Yxsof   = Par(16);  %biomass yield from the overflow route
 % g o2 L-1 monod O0
% qOmax = 13.4*32/1000;
qS      = (qSmax/(1+(A/Kia)))*(S/(S+Ks)); % Check
qSof    = pAmax*(qS/(qS+Kap));             % Check
pA      = qSof*Yas;                        % Check
qSox    = (qS-qSof)*(O/(O+Ko));        % Check      
%qSan    = (qSox-qm)*Yem*Cx/Cs;             % 
qsA     = (qAmax/(1+(qS/Kis)))*(A/(A+Ksa));
qA      = pA-qsA;                            

qO      = Yos*(qSox-qm)+qsA*Yoa;

mu      = (qSox-qm)*Yem + qsA*Yxa + qSof*Yxsof;

sys = [X S A O V mu];
sys(sys < 0) = 0;
end