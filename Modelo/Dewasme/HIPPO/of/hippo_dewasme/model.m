
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
dxdt = [];
% Model parameters:
% State variables
X = x(1);   % Biomass (g/L)
S = x(2);   % Substrate (g/L)
E = x(3);   % Ethanol (g/L)
O = x(4);   % Dissolved Oxygen (g/L)


% Fixed parameters
muset   = 0.197;
X0      = 4.04;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.035;
Ko      = 0.0001;
V  = V0;

% Adjusted parameters (all cases)
% Ks      = k(1);
qSmax   = k(1);
Ysoxx   = k(2);
% Yso     = k(4);
% Kio     = k(5);
Yse     = k(3);
% Kec     = k(7);
% Ysofx   = k(8);
% Yeo     = k(9);
% Yex     = k(10);
% qOmax   = k(6);
% Yosof   = k(12);
Ysofx = k(4); 
% Yse = k(8);

% FIXED
Kec     = 0.0051;
Yeo     = 1.998;
Yex     = 1.7531;
qOmax   = 0.3497;
Yosof   = 0.0001;
Yso     = 0.3438;
Kio     = 4.7323;
Ks      = 0.1414;


% % Adjusted parameters (overflow)
% Kie     = kof(1);
% Yes     = kof(2);
% Kec     = kof(3);
% Ysofx   = kof(4);
% Yoe     = kof(5);
% Yxe     = kof(6);
% qOmax   = kof(7);
% Yosof   = kof(8);
a = 100;
% Constitutive equations
qS      = qSmax*S/(S+Ks);
qO      = (qOmax*O/(O+Ko))*(Kio/(Kio+E));
qScrit  = qO/Yso;
qSox    = min(qS,qScrit);
% qSox    = (qS*exp(-a*qS)+qScrit*exp(-a*qScrit))/(exp(-a*qS)+exp(-a*qScrit));
qSof    = max(0,qS-qScrit);
% qSof    = (0*exp(a*0)+(qS-qScrit)*exp(a*(qS-qScrit)))/(exp(a*0)+exp(a*(qS-qScrit)));
% qec     = (Yso/Yeo)*(qScrit-qS)*(E/(E+Kec));
% qE      = (0*exp(a*0)+qec*exp(a*qec))/(exp(a*0)+exp(a*qec));
qE      = max(0,(Yso/Yeo)*(qScrit-qS)*(E/(E+Kec)));
mu      = Ysoxx*qSox + Ysofx*qSof + Yex*qE;
F0      = muset*(X0*V0)/(Ysoxx*Sin);
F       = F0*exp(muset*t);
D       = F/V;

% ODEs
dxdt(1) = X*(mu-D);                 % dXdt
dxdt(2) = D*(Sin-S)-(qSox-qSof)*X;  % dSdt
dxdt(3) = (Yse*qSof-qE)*X-D*E;      % dEdt
dxdt(4) = klao2*(osat-O)-D*O - (Yso*qSox+Yosof*qSof+Yeo*qE)*X; % dOdt
                 

dxdt = dxdt';



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%