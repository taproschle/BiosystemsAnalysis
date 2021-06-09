function dydt = dewasme_unified(t,y,v,k,kof)

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
Ks      = k(1);
qSmax   = k(2);
Ysoxx   = k(3);
Yos     = k(4);

% Adjusted parameters (overflow)
Kie     = kof(1);
Yes     = kof(2);
Kec     = kof(3);
Ysofx   = kof(4);
Yoe     = kof(5);
Yxe     = kof(6);
qOmax   = kof(7);
Yosof   = kof(8);

% Constitutive equations
qS      = qSmax*S/(S+Ks);
qO      = (qOmax*O/(O+Ko))*(Kie/(Kie+E));
qScrit  = qO/Yos;
qSox    = min(qS,qScrit);
qSof    = max(0,qS-qScrit);
qE      = max(0,(Yos/Yoe)*(qScrit-qS)*(E/(E+Kec)));
mu      = Ysoxx*qSox + Ysofx*qSof + Yxe*qE;
F0      = muset*(X0*V0)/(Ysoxx*Sin);
F       = F0*exp(muset*t);
D       = F/V;

% ODEs
dydt(1) = X*(mu-D);                 % dXdt
dydt(2) = D*(Sin-S)-(qSox-qSof)*X;  % dSdt
dydt(3) = (Yes*qSof-qE)*X-D*E;      % dEdt
dydt(4) = klao2*(osat-O)-D*O-...    % dOdt
    (Yos*qSox+Yosof*qSof+Yoe*qE)*X;
dydt(5) = F;                        % dVdt

dydt = dydt';

end