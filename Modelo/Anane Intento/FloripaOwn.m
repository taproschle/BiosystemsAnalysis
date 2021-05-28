function [sys, x0,str,ts] = FloripaOwn(t,x,u,flag)

switch flag
   case 0
   [sys,x0,str,ts] = mdlInitializeSizes;
   case 1 
   sys = mdlDerivatives(t,x,u);
   case 3
   sys = mdlOutputs(t,x,u);
   case {2,4,9}
   sys=[];
   otherwise
   error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts] = mdlInitializeSizes
sizes=simsizes;
sizes.NumContStates  = 6;  % (biomass, glucose, ethanol, O2, CO2, volume) 
sizes.NumDiscStates  = 0;  
sizes.NumOutputs     = 9;  % (biomass, glucose, ethanol, O2, CO2, volume)
sizes.NumInputs      = 4;  % (substrate flowrate, agitation speed, gas flowrate)
sizes.DirFeedthrough = 0; 
sizes.NumSampleTimes = 1; 
sys = simsizes(sizes);
str = []; 
ts = [0 0];
%load 'CondicionesIniciales.mat '  C_iniciales
x0 = C_iniciales;     % Initial conditions

function sys = mdlDerivatives(~,x,u)
global  kx1 kx2 kx3 ks1 ks2 ke2 ke3 ko1 ko2 ko3 kc1 kc2 kc3 kos koe ...
        rs_max muo Ko  Ks  Ke  Kie  CO2s alfa beta gamma delta Henry Si P
% State variables
X   =  x(1);                     % Biomas,           [g/L]
S   =  x(2);                     % Glucose,          [g/L]
E   =  x(3);                     % Ethanol,          [g/L]
O2  =  x(4);                     % Dissolved O2,     [g/L]
CO2 =  x(5);                     % Dissolved CO2,    [g/L]
V   =  x(6);                     % Volume,           [L]

Fs  = u(1);                      % Feed rate,        [L/h]
N   = u(2);                      % Agitation rate    [rpm]
G   = u(3);                      % Aire flow         [LPM]
yo2 = u(4);                      % Fraccion gaseosa de O2

% =========================== Algebraic equations =========

% Mass transfer coefficient
klao2 = 5*alfa*(N^beta)*(G^gamma);                   % Mass transfer coefficient for O2,[1/h]
klaco2 = delta*klao2;   


% ============================ Restrictions 
%                                              
r1 = min(rsu,rscrit)/ks1;                            % Glucose oxidation rate,       [g of S/g of X/h]
r2 = max(0,rsu-rscrit)/ks2;                          % Glucose fermentation rate     [g of S/g of X/h]
r3 = max(0, kos*(rscrit-rsu)*(E/(E+Ke))/koe)/ke3;    % Ethanol oxidation rate        [g of E/g of X/h]  
 
mu = (kx1*r1 + kx2*r2 + kx3*r3);                     % Specific growth rate,                           [1/h]
rs = (ks1*r1 + ks2*r2);                              % Specific glucose consumption rate [1/h]
ra = (ke2*r2 - ke3*r3);                              % Specific ethanol rate,                          [1/h]
ro = (ko1*r1 + ko2*r2 + ko3*r3);                     % Specific O2 consumption rate,     [1/h]
rc = (kc1*r1 + kc2*r2 + kc3*r3);                     % Specific CO2 consumption rate,    [1/h]

% ========================== Differential equations
sys(1) = X*mu - X*D;                                 % Biomass balance,        [g/L/h]
sys(2) = Si*D - S*D - rs*X;                          % Glucose balance,        [g/L/h]  
sys(3) = X*ra - E*D;                                 % Ethanol balance,        [g/L/h]
sys(4) = klao2*(P*yo2/Henry-O2) - ro*X - O2*D;       % O2 balance,             [g/L/h]
sys(5) = -klaco2*(CO2-CO2s) + rc*X - CO2*D;          % CO2 balance,            [g/L/h] 
sys(6) = Fs;                                         % Total mass balance,     [L/h]


function sys = mdlOutputs(~,x,~)
global  kx1 kx2 kx3 ks1 ks2 ke3 kos koe rs_max muo Ko  Ks  Ke  Kie

% State variables
%X   =  x(1);                      % Biomas,           [g/L]
S   =  x(2);                       % Glucose,          [g/L]
E   =  x(3);                       % Ethanol,          [g/L]
O2  =  x(4);                       % Dissolved O2,     [g/L]
%CO2 =  x(5);                      % Dissolved CO2,    [g/L]
%V   =  x(6);                      % Volume,           [L]

% Algebraic equations
rsu = rs_max*(S/(Ks+S));                              % Glucose unlimited rate, [g of S/g of X/h]
rscrit = muo*(O2/(Ko+O2))*(Kie/(Kie+E))/kos;          % Critical glucose oxidation  rate [g of S/g of X/h]
% Restrictions
r1 = min(rsu,rscrit)/ks1;                             % Glucose oxidation rate,    [g of S/g of X/h]
r2 = max(0,rsu-rscrit)/ks2;                           % Glucose fermentation rate, [g of S/g of X/h]
r3 = max(0, kos*(rscrit-rsu)*(E/(E+Ke))/koe)/ke3;     % Ethanol oxidation rate [g of E/g of X/h]  
mu   = (kx1*r1 + kx2*r2 + kx3*r3);                    % Specific growth rate [1/h]
Scrit = Ks * rscrit / ( rs_max - rscrit);             % Critical substrate [g/L]

sys = [x(1) x(2) x(3) x(4) x(5) x(6) mu rscrit*kx1 Scrit]; 









