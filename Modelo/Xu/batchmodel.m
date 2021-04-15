function dydt = batchmodel(t,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dydt = [];

% parametros
Yas = 0.23; %  g A /g S stochiometric cte CARCAMO
%Yas = 0.667; % XU
Yoa = 20.66/1000; % g o2 / g A CARCAMO
%Yoa = 1.067; % XU 

Yos = 401/1000; % g o2 / g S CARCAMO
%Yos = 1.067; % XU
%Yxa = 14.55; % g DCW / gIP CARCAMO
Yxa = 0.4; % XU
Yxsof = 0.7; %CARCAMO
%Yxsof = 0.15; % XU
Yxsox = 0.3; % CARCAMO
%Yxsox = 0.51; %XU


%%
Ca = 1/30; % mol Carbon / g acetato
Cs = 1/30; % mol Carbon /g azucar
Cx = 0.04; % mol Carbon / g biomasa


%%
%KiO = 4; % g/L fitting
Ka = 0.5236; % CARCAMO
%Ka = 0.05; %XU
KiS = 0.45; % g A /L fitting CARCAMO
%KiS = 5; % XU
Ks = 8; % g S/L CARCAMO
%Ks = 0.05; % XU
Ko = 0.0045; %g O/L CARCAMO
%Ko = 1/100;
qAcmax = 6; % g/g hr CARCAMO
%qAcmax = 0.06; % XU
qm = 0.04; % g/ g hr
qOmax = 83.5/1000; % g/g hr CARCAMO
%qOmax = 15.6*32/1000; % XU
qSmax = 4.9; % g/g hr CARCAMO
%qSmax = 1.3; % Xu
muc = 0.55;
Osat = 0.035;
kla = 180; % [h-1]


% variables de estado
S = x(1); % sustrato [ g/L]
%disp(S)
A = x(2); % acetato [g/L]
X = x(3); % biomasa [g/L]
V = x(4); % volumen hasta ahora son 100ml parece [L]
O = x(5); %g/L
% ecs
% oxidative pathway; % glucose uptake [g g-1 celulas h-1]
% relacion rutas metabolicas
qS_T = qSmax*S / ((Ks + S)*(1 + A/KiS));
qO_crit = qOmax*O / (Ko + O);
qS_ox = min (qS_T , qO_crit /Yos ) ; 
qS_of   = max (0 , qS_T - qS_ox) ; 
qAc = min (qAcmax * A / (A + Ka) , (qO_crit - qS_ox*Yos) / Yoa) ;
% ecs constitutivas 
qS = qS_ox + qS_of;
qA = qS_of*Yas - qAc; % qA produccion - qA de consumo
mu = qS_ox*Yxsox + qS_of*Yxsof + qAc*Yxa; 
qO2 = qS_ox*Yos + qAc*Yoa; % lo q se consume en la via oxidativa + lo que se consume en la via overflow

%qO = qOs + qAc*Yoa;

%mu = (qS_ox - qm)*Yxsox + qS_of*Yxsof + qAc*Yxa;

%edo

dydt(1) =  - qS*X; %sustrato
%if t>=12 && t<=12.5
    %dydt(1) = -qS*X + 10;
%end
dydt(2) = (qA)*X; % acetato
%dydt(2) = (mu)*X;
dydt(3) = (mu)*X; % biomass
dydt(4) = 0; % volumen
dydt(5) = -(qO2)*X + kla*(Osat-O);
dydt = dydt';



end

