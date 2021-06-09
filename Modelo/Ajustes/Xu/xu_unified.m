function dydt = xu_unified(t,y,v,k,kof)

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
Yos     = k(5);

% Adjusted parameters (overflow)
Kie     = kof(1);
Yes     = kof(2);
Kec     = kof(3);
qEmax   = kof(4);
Ysofx   = kof(5);
Yoe     = kof(6);
Yxe     = kof(7);
qOmax   = kof(8);

% Constitutive equations
qS      = (qSmax*S/(Ks+S))*1/(1+(E/Kie));
qOs     = min(qOmax*(O/(Ko+O))*(1/(1+(E/Kie))),qOmax);
qSox    = min(((qOs/Yos)-(qm*Ysoxx))/(1-Ysoxx),qS);
qSof    = max(qS-qSox,0);
qEp     = (qSof-qSof*Ysofx)*Yes;
qEc     = min(qEmax*E/(E+Kec),(qOmax-qOs)*Yoe/(1-Yxe));
mu      = (qSox-qm)*Ysoxx+qSof*Ysofx+qEc*Yxe;
qO      = qOs+(qEc-qEc*Yxe)*Yoe;
F0      = muset*(X0*V0)/(Ysoxx*Sin);
F       = F0*exp(muset*t);
D       = F/V;

% ODEs
dydt(1) = X*(mu-D);                 % dXdt
dydt(2) = D*(Sin-S)-qS*X;           % dSdt
dydt(3) = (qEp-qEc)*X-D*E;          % dEdt
dydt(4) = klao2*(osat-O)-qO*X-D*O;  % dOdt
dydt(5) = F;                        % dVdt

dydt = dydt';

end