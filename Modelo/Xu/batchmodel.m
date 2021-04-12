function dydt = batchmodel(t,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dydt = [];

% parametros
Yas = 0.667; % g/g stochiometric cte
Yoa = 1.067;
Yos = 1.067;
Yxa = 0.4;
Yxsof = 0.15;
Yxsox = 0.51;
Ca = 1/30; % mol Carbon / g acetato
Cs = 1/30; % mol Carbon /g azucar
Cx = 0.04; % mol Carbon / g biomasa
Ka = 0.05; % g/L
KiO = 4; % g/L fitting
KiS = 5; % g/L fitting
Ks = 0.05; % g/L
qAcmax = 0.06; % g/g hr
qm = 0.04; % g/ g hr
qOmax = 15.6; % mmol/g hr
qSmax = 1.3; % g/g hr
muc = 0.55;

% variables de estado
S = x(1); % sustrato [ g/L]
%disp(S)
A = x(2); % acetato [g/L]
X = x(3); % biomasa [g/L]
V = x(4); % volumen hasta ahora son 100ml parece [L]
% ecs
% oxidative pathway
qS = ((qSmax)/(1+(A/KiS)))*(S/(S+Ks)); % glucose uptake [g g-1 celulas h-1]
%disp(qS)
qSox = qS;
qSox_an = (qSox - qm)*Yxsox*Cx/Cs; 
qSox_en = qSox - qSox_an;
qOs = qSox_en*Yos;
%qSox = qSox_an + qSox_en;
%qSof = qS - qSox; 
%disp(qOs)
%if qOs > qOmax/(1+(A/KiO))-100
    %qSof = qS*0.15;
%end
%qSof = qS*0.25;
qSof = qS - qSox;
%overflow
qSof_an = qSof*Yxsof*Cx/Cs;
qSof_en = qSof - qSof_an;
qAp = qSof_en*Yas;
qAc = qAcmax*A/(A + Ka);
%disp(qAc)
qAc_an = qAc*Yxa*Cx/Ca;
qAc_en = qAc - qAc_an;
%if qAc_en > ((qOmax - qOs)/Yoa)
    %qAc_en = ((qOmax - qOs)/Yoa);
%end

qO = qOs + qAc_en*Yoa;
%disp(qO)
%if t < 8
    %qO = 12;
%end


%disp(qO);
mu = (qSox - qm)*Yxsox + qSof*Yxsof + qAc*Yxa;

%edo

dydt(1) =  - qS*X; %sustrato
%if t>=12 && t<=12.5
    %dydt(1) = -qS*X + 10;
%end
%dydt(2) = (qAp - qAc)*X; % acetato
dydt(2) = (mu-muc)*X;
dydt(3) = (muc)*X; % biomass
dydt(4) = 0; % volumen

dydt = dydt';



end

