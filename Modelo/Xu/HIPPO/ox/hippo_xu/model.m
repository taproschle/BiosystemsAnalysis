
function dxdt = model(t,x,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DO NOT MODIFY THIS SECTION

kfixed = evalin('base','kfixed');
p = zeros(length(kfixed));

j = 1;
for i = 1:length(kfixed)
    if isnan(kfixed(i))
        p(i) = k(j);
        j    = j+1;
    else
        p(i) = kfixed(i);
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EDIT HERE

% Model parameters:
% State variables
X = x(1);   % Biomass (g/L)
S = x(2);   % Substrate (g/L)
E = x(3);   % Ethanol (g/L)
O = x(4);   % Dissolved Oxygen (g/L)


% Fixed parameters
muset   = 0.11;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.008;
Ko      = 0.0001;
Kio     = 4;
V = V0;

% Adjusted parameters (all cases)
% Ks      = k(5);
qSmax   = k(1);
Ysoxx   = k(3);
% qm      = k(4);
Yso     = k(5);
% Kie     = k(7);
% Yse     = k(13);
% Kec     = k(6);
% qEcmax   = k(2);
Ysofx   = k(4);
% Yeo     = k(11);
% Yex     = k(12);
qOmax   = k(2);

Yeo     = 0.7607;
Yex     = 0.155;
Kie     = 7.4378;
Yse     = 0.3652;
Kec     = 0.5019;
qEcmax   = 0.1224;
qm      = 0.0575;
Ks      = 0.1034;

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
a = 100;
% Constitutive equations
qS      = (qSmax*S/(Ks+S))*1/(1+(E/Kie));
qOs     = qOmax*(O/(Ko+O))*(1/(1+(E/Kio)));
qScrit = ((qOs/Yso)-(qm*Ysoxx))/(1-Ysoxx);
% qSox    = (qS*exp(-a*qS)+qScrit*exp(-a*qScrit))/(exp(-a*qS)+exp(-a*qScrit));
% qSof    = (0*exp(a*0)+(qS-qSox)*exp(a*(qS-qSox)))/(exp(a*0)+exp(a*(qS-qSox)));
% qeci = qEcmax*E/(E+Kec);
% qecj = (qOmax-qOs)*Yeo/(1-Yex);
qSox    = min(qS,qScrit);
qSof    = max(qS-qSox,0);
qEp     = (qSof-qSof*Ysofx)*Yse;
% qEc      = (qeci*exp(-a*qeci)+qecj*exp(-a*qecj))/(exp(-a*qeci)+exp(-a*qecj));
qEc     = min(qEcmax*E/(E+Kec),(qOmax-qOs)*Yeo/(1-Yex));
mu      = (qSox-qm)*Ysoxx+qSof*Ysofx+qEc*Yex;
qO      = qOs+(qEc-qEc*Yex)*Yeo;
F0      = muset*(X0*V0)/(Ysoxx*Sin);
F       = F0*exp(muset*t);
D       = F/V;

% ODEs
dxdt(1) = X*(mu-D);                 % dXdt
dxdt(2) = D*(Sin-S)-qS*X;           % dSdt
dxdt(3) = (qEp-qEc)*X-D*E;          % dEdt
dxdt(4) = klao2*(osat-O)-qO*X-D*O;  % dOdt
                      % dVdt

dxdt = dxdt';


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%