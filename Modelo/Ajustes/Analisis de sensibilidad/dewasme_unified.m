
function dydt = dewasme_unified(t,y,k)

% State variables
X = y(1);   % Biomass (g/L)
S = y(2);   % Substrate (g/L)
E = y(3);   % Ethanol (g/L)
O = y(4);   % Dissolved Oxygen (g/L)
V = y(5);   % Volume (L)

% Fixed parameters
muset   = 0.11;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.035;
Ko      = 0.0001;


% Adjusted parameters (all cases)
Ks      = k(1);
qSmax   = k(2);
Ysoxx   = k(3);
Yso     = k(4);
Kio     = k(5);
Yse     = k(6);
Kec     = k(7);
Ysofx   = k(8);
Yeo     = k(9);
Yex     = k(10);
qOmax   = k(11);
Yosof   = k(12);

% % Adjusted parameters (overflow)
% Kie     = kof(1);
% Yes     = kof(2);
% Kec     = kof(3);
% Ysofx   = kof(4);
% Yoe     = kof(5);
% Yxe     = kof(6);
% qOmax   = kof(7);
% Yosof   = kof(8);

% Constitutive equations
qS      = qSmax*S/(S+Ks);
qO      = (qOmax*O/(O+Ko))*(Kio/(Kio+E));
qScrit  = qO/Yso;
qSox    = min(qS,qScrit);
qSof    = max(0,qS-qScrit);
qE      = max(0,(Yso/Yeo)*(qScrit-qS)*(E/(E+Kec)));
mu      = Ysoxx*qSox + Ysofx*qSof + Yex*qE;
F0      = muset*(X0*V0)/(Ysoxx*Sin);
F       = F0*exp(muset*t);
D       = F/V;

% ODEs
dydt(1) = X*(mu-D);                 % dXdt
dydt(2) = D*(Sin-S)-(qSox-qSof)*X;  % dSdt
dydt(3) = (Yse*qSof-qE)*X-D*E;      % dEdt
dydt(4) = klao2*(osat-O)-D*O-...    % dOdt
    (Yso*qSox+Yosof*qSof+Yeo*qE)*X;
dydt(5) = F;                        % dVdt

dydt = dydt';

end