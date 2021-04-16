function dY = modelo_carcamo(t,x)
global Vin Xin Factor_mu YS_ox_X YS_ox_O YS_f_X YS_f_IP mu_set YIP_X YIP_O qS_max qIP_max qO_max K_S K_IP Ki_IP K_O kLa_O2 Sfeed DO_max mu_max
% Parametros
% Coeficientes de rendimiento
% Regimen oxidativo - altamente energético  32 ATP
% YS_ox_X = p(1); % [gDCW / g_Glu] biomasa en el consumo oxidativo de glucosa
% YS_ox_O = p(2); % [mg_O2 / g_Glu] oxígeno en el consumo oxidativo de glucosa
% % Regimen fermentativo - 2 ATP
% YS_f_X = p(3); % [gDCW / g_Glu] biomasa en el consumo fermentativo de glucosa
% YS_f_IP = p(4); % [g_IP / g_Glu ]
% % Reaccion de Subproducto
% YIP_X = p(5); % [g_DCW / g_IP ] biomasa en el consumo oxidativo del producto inhibitorio
% YIP_O = p(6); % [mg_O2 / g_IP ] oxígeno en el consumo oxidativo del producto inhibitorio
% % Velocidades específicas de reaccion
% qS_max = p(7); % [g_Glu / (g_DCW * h)] tasa específica de consumo máxima de glucosa
% qIP_max = p(8); % [g_Glu / (g_IP * h)]
% qO_max = p(9); % [mg_O2 / (g_DCW * h)] tasa específica de consumo de oxígeno máxima
 qS_crit = qO_max / YS_ox_O ; %  [g_Glu / (g_DCW * h)]
% tasa específica de crecimiento crítica
mu_crit = qS_crit * YS_ox_X ; % [1/h]

kLa_O2 = 140 * 20 ; % [1/h] 
DO_max = 35; % [mg_O2 / L] Ley de Henry *** Concentración de oxígeno maxima
% Sin = 500; % Sustrato en la corriente de alimentacion - Parametro del balance de sustrato
%mu_set = mu_crit * Factor_mu;

% Asignacion de variables
X = x(1);
S = x(2);
P = x(3);
V = x(4);
DO = x(5);

%F = mu_set*Vin*Xin*exp(mu_set*t) / (YS_ox_X*Sfeed);

F = mu_set*Vin*Xin*exp(mu_set*t )/(YS_ox_X*Sfeed);
%disp(t)


% Relaciones constitutivas:
q_ST = qS_max*S / ((K_S + S)*(1 + P/Ki_IP));
qO_crit = qO_max*DO / (K_O + DO);

qS_ox = min (q_ST , qO_crit / YS_ox_O ) ; 
qS_f   = max (0 , q_ST - qS_ox) ; 
qIP_cons = min (qIP_max * P / (P + K_IP) , (qO_crit - qS_ox*YS_ox_O) / YIP_O) ;

% Ecuaciones constitutivas 

qS = qS_ox + qS_f;
qIP = qS_f*YS_f_IP - qIP_cons;
%mu = qS_ox*YS_ox_X + qS_f*YS_f_X + qIP_cons*YIP_X; 
 mu = min(qS_ox*YS_ox_X + qS_f*YS_f_X + qIP_cons*YIP_X,mu_max);

qO2 = qS_ox*YS_ox_O + qIP_cons*YIP_O;

% Ecuaciones diferenciales:
dY = zeros (5,1);
% Biomasa
% d(XV) / dt
dY (1) = mu*X -(F/V)*X;

% Sustrato
% d(SV) / dt
dY (2) = -qS*X + F*(Sfeed)/V;

% Producto
% d(IP V) / dt
dY (3) =  qIP*X;
% Volumen
% dV / dt
dY (4) = F;

% Oxigeno
% d(DO V) / dt
dY (5) = kLa_O2 * (DO_max - DO) - qO2*X ;