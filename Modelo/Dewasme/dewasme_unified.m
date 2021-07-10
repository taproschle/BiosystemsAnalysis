
function dydt = dewasme_unified(t,y,v,k)

% State variables
X = y(1);   % Biomass (g/L)
S = y(2);   % Substrate (g/L)
E = y(3);   % Ethanol (g/L)
O = y(4);   % Dissolved Oxygen (g/L)
V = y(5);   % Volume (L)

% Fixed parameters
muset   = v(1);
X0      = v(2);
V0      = v(3);
Sin     = v(4);
klao2   = v(5);
osat    = v(6);
Ko      = v(7);

% Adjusted parameters (all cases)
qSmax = k(1);
Ysoxx = k(2);
Kec = k(3);
Yex = k(4);
qOmax = k(5);


%FIXED
% qSmax = 1;
Yeo     = 1.998;
Yosof   = 0.0001;
Yso     = 0.3438;
Kio     = 4.7323;
Ks      = 0.1414;
Yse     = 0.9634;
Ysofx   = 0.3294; 
% Ysofx = 0;

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