function dydt = xu_model(t,y,v,k)

% PARAMETROS AJUSTABLES
YS_ox_X = k(1); YS_of_X = k(2); Yax = k(3); 
Yoa = k(4); Ysa = k(5); Yso = k(6); 
C_X = k(7) ; C_A = k(8); C_S = k(9);
qm = k(10); qS_max = k(11); qAc_max = k(12); qO_max = k(13);  
K_S = k(14); K_i_A = k(15); K_A = (16);
 
% PARAMETROS FIJOS
mu_set  = v(1); 
klao2   = v(2); 
Vin     = v(3);
Xin     = v(4);
Sfeed   = v(5);
O_sat   = v(6);
K_O     = v(7);

dydt = [];

% variables de estado
X = y(1); % biomasa [ g/L]
S = y(2); % glucosa [g/L]
A = y(3); % acetato [g/L]
O = y(4); % oxigeno disuelto [g/L]
V = y(5); % volumen

%ec constitutivas
qS = (qS_max*S/(K_S + S))*1/(1+(A/K_i_A));
qOs = min(qO_max*(O/(K_O+O))*(1/(1+(A/K_i_A))),qO_max);
qSox = min(((qOs/Yso)-(qm*YS_ox_X*C_X/C_S))/(1-YS_ox_X*C_X/C_S),qS);
qSof = max(qS - qSox,0);
qAp = (qSof - qSof*YS_of_X*C_X/C_S)*Ysa;
qAc = min(qAc_max*A/(A+K_A) , (qO_max - qOs)*Yoa/(1-Yax*C_X/C_A));

mu = (qSox - qm)*YS_ox_X + qSof*YS_of_X + qAc*Yax;
qO = qOs + (qAc - qAc*Yax*C_X/C_A)*Yoa;

F = mu_set*Vin*Xin*exp(mu_set*t)/(YS_ox_X*Sfeed);

% Ecuaciones diferenciales.
dydt(1) = (mu)*X - (F/V)*X; % biomass
dydt(2) =  (F/V)*(Sfeed-S)- qS*X ; %sustrato
%if t>=12 && t<=12.5
    %dydt(1) = -qS*X + 10;
%end
dydt(3) = (qAp-qAc)*X - (F/V)*A; % acetato
%dydt(2) = (mu)*X;
dydt(4) = klao2*(O_sat - O)-qO*X - (F/V)*O;
dydt(5) = F; % volumen

dydt = dydt';
end
