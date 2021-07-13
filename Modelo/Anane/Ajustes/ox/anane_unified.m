function dydt = anane_unified(t,y,v,k)

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

% Fixed parameters (all cases)
Ks      = 0.0450;
qm      = 0.0067;
Yso     = 0.4684;
Kie     = 1.6129;
qEpmax  = 0.6350;
Kep     = 0.5034;
Yse     = 0.1955;
Kis     = 2.2729;
Yeo     = 0.8163;

% Adjusted parameters (all cases)
qSmax   = k(1);
qEcmax  = k(2);
Kec     = k(3);
Ysoxx   = k(4);
Ysofx   = k(5);
Yex     = k(6);


% Constitutive equations
qS      = (qSmax/(1+(E/Kie)))*(S/(S+Ks));
qSof    = qEpmax*(qS/(qS+Kep));
qEp      = qSof*Yse;
qSox    = (qS-qSof)*(O/(O+Ko));
qEc     = (qEcmax/(1+(qS/Kis)))*(E/(E+Kec));
qE      = qEp-qEc;
mu      = (qSox-qm)*Ysoxx+qSof*Ysofx+qEc*Yex;
qO      = Yso*(qSox-qm)+qEc*Yeo;
F0      = muset*(X0*V0)/(Ysoxx*Sin);
F       = F0*exp(muset*t);
D       = F/V;

% ODEs
dydt(1) = X*(mu-D);                 % dXdt
dydt(2) = D*(Sin-S)-qS*X;           % dSdt
dydt(3) = qE*X-D*E;                 % dEdt
dydt(4) = klao2*(osat-O)-qO*X-D*O;  % dOdt
dydt(5) = F;                        % dVdt

dydt = dydt';

end