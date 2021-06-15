function dydt = xu_unified(t,y,k)

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
osat    = 0.008;
Ko      = 0.0001;
Kio     = 4;


% Adjusted parameters (all cases)
Ks      = k(1);
qSmax   = k(2);
Ysoxx   = k(3);
qm      = k(4);
Yso     = k(5);
Kie     = k(6);
Yse     = k(7);
Kec     = k(8);
qEcmax   = k(9);
Ysofx   = k(10);
Yeo     = k(11);
Yex     = k(12);
qOmax   = k(13);


% % Adjusted parameters (overflow)
% Kie     = kof(1);
% Yes     = kof(2);
% Kec     = kof(3);
% qEmax   = kof(4);
% Ysofx   = kof(5);
% Yoe     = kof(6);
% Yxe     = kof(7);
% qOmax   = kof(8);
% Kio     = kof(9);

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
dydt(5) = F;                        % dVdt

dydt = dydt';

end