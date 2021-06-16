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

% Adjusted parameters (all cases)
Ks      = k(1);
qSmax   = k(2);
Ysoxx   = k(3);
qm      = k(4);
Yso     = k(5);
Kie     = k(6);
qEpmax   = k(7);
Kep     = k(8);
Yse     = k(9);
Kec     = k(10);
qEcmax   = k(11);
Kis     = k(12);
Ysofx   = k(13);
Yeo     = k(14);
Yex     = k(15);

% % Adjusted parameters (overflow)
% Kie     = kof(1);
% pEmax   = kof(2);
% Kep     = kof(3);
% Yes     = kof(4);
% Kec     = kof(5);
% qEmax   = kof(6);
% Kis     = kof(6);
% Ysofx   = kof(7);
% Yoe     = kof(8);
% Yxe     = kof(9);

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