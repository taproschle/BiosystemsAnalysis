function dydt = fedbatch(t,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global klao2 O_sat K_O Yax YS_ox_X YS_of_X Yoa K_A C_X C_A C_S qm qO_max qAc_max qS_max mu_set Ysa Yso K_i_A K_S Vin Xin Sfeed
dydt = [];

% variables de estado
S = x(1); % sustrato [ g/L]
%disp(S)
A = x(2); % acetato [g/L]
X = x(3); % biomasa [g/L]
V = x(4); % volumen  hasta ah
O = x(5); % oxigeno disuelto [g/L]
%ec constitutivas

qS = (qS_max*S/(K_S + S))*1/(1+(A/K_i_A));
%qOs = min((qS-(qS-qm)*YS_ox_X*C_X/C_S)*Yso , qO_max); 

%qOs = min(qO_max*(O/(K_O+O))*(1/(1+(A/K_i_A))),qO_max*(1/(1+(A/K_i_A))));
qOs = min(qO_max*(O/(K_O+O))*(1/(1+(A/K_i_A))),qO_max);
qSox = min(((qOs/Yso)-(qm*YS_ox_X*C_X/C_S))/(1-YS_ox_X*C_X/C_S),qS);
%qSox = min(qOs/Yso,qS);
%qSox = qOs/Yso;
qSof = max(qS - qSox,0);
qAp = (qSof - qSof*YS_of_X*C_X/C_S)*Ysa;
qAc = min(qAc_max*A/(A+K_A) , (qO_max - qOs)*Yoa/(1-Yax*C_X/C_A));

mu = (qSox - qm)*YS_ox_X + qSof*YS_of_X + qAc*Yax;
qO = qOs + (qAc - qAc*Yax*C_X/C_A)*Yoa;
    
F = mu_set*Vin*Xin*exp(mu_set*t)/(YS_ox_X*Sfeed);



dydt(1) =  (F/V)*(Sfeed-S)- qS*X ; %sustrato
%if t>=12 && t<=12.5
    %dydt(1) = -qS*X + 10;
%end
dydt(2) = (qAp-qAc)*X - (F/V)*A; % acetato
%dydt(2) = (mu)*X;
dydt(3) = (mu)*X - (F/V)*X; % biomass
dydt(4) = F; % volumen
dydt(5) = klao2*(O_sat - O)-qO*X - (F/V)*O;
dydt = dydt';

end
