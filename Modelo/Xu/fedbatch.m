function dydt = fedbatch(t,x,u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global Yas Yoa Yos Yxa Yxsof Yxsox Ka KiS Ks Ko qAcmax qOmax qSmax Osat kla mumax
dydt = [];
mu_set = u(1);
Xin = u(2);
Vin = u(3);
Sin = u(4);

% parametros



%%
Ca = 1/30; % mol Carbon / g acetato
Cs = 1/30; % mol Carbon /g azucar
Cx = 0.04; % mol Carbon / g biomasa


%%


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
% if mu >= mumax
%     mu = mumax;
% end
% disp(mu)
qO2 = qS_ox*Yos + qAc*Yoa; % lo q se consume en la via oxidativa + lo que se consume en la via overflow

%qO = qOs + qAc*Yoa;

%mu = (qS_ox - qm)*Yxsox + qS_of*Yxsof + qAc*Yxa;

    
F = mu_set*Vin*Xin*exp(mu_set*t)/(Yxsox*Sin);


dydt(1) =  (F/V)*(Sin-S)- qS*X ; %sustrato
%if t>=12 && t<=12.5
    %dydt(1) = -qS*X + 10;
%end
dydt(2) = (qA)*X - (F/V)*A; % acetato
%dydt(2) = (mu)*X;
dydt(3) = (mu)*X - (F/V)*X; % biomass
dydt(4) = F; % volumen
dydt(5) = -(qO2)*X + kla*(Osat-O);
dydt = dydt';

end
