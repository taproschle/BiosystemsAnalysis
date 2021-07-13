function dydt = xu_unified(t,y,v,k)

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
Kio     = v(8);

% Adjusted parameters (all cases) %0.239492        0.2118       0.730404         0.846497
% qSmax   = k(1);
qSmax   = k(1);
qOmax   = k(2); 
Ysoxx   = k(3);
Ysofx   = k(4);
Yso     = k(5);
% FIX
% Yso = 0.8965;
Yeo     = 0.7607;
Yex     = 0.155;
Kie     = 7.4378;
Yse     = 0.3652;
Kec     = 0.5019;
qEcmax   = 0.1224;
qm      = 0.0575;
Ks      = 0.1034;


% Constitutive equations
qS      = (qSmax*S/(Ks+S))*1/(1+(E/Kie));
qOs     = min(qOmax*(O/(Ko+O))*(1/(1+(E/Kio))),qOmax);
qSox    = min(((qOs/Yso)-(qm*Ysoxx))/(1-Ysoxx),qS);
qSof    = max(qS-qSox,0);
qEp     = (qSof-qSof*Ysofx)*Yse;
qEc     = min(qEcmax*E/(E+Kec),(qOmax-qOs)*Yeo/(1-Yex));
mu      = (qSox-qm)*Ysoxx+qSof*Ysofx+qEc*Yex;
qO      = qOs+(qEc-qEc*Yex)*Yeo;
F0      = muset*(X0*V0)/(Ysoxx*Sin);
F       = F0*exp(muset*t);
D       = F/V;

% ODEs
dydt(1) = X*(mu-D);                 % dXdt
dydt(2) = D*(Sin-S)-qS*X;           % dSdt
dydt(3) = (qEp-qEc)*X-D*E;          % dEdt
dydt(4) = klao2*(osat-O)-qO*X-D*O;  % dOdt
% dydt(4) = 0;
dydt(5) = F;                        % dVdt

dydt = dydt';

end