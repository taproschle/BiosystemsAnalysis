
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
V = x(5);

% Fixed parameters
muset   = 0.11;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.008;
Ko      = 0.0001;


% Adjusted parameters (all cases)
% Ks      = k(5);
qSmax   = p(1);
Ysoxx   = p(4);
% qm      = k(4);
% Yso     = k(12);
% Kie     = k(8);
% qEpmax   = k(3);
% Kep     = k(8);
% Yse     = k(15);
Kec     = p(3);
qEcmax   = p(2);
% Kis     = k(9);
Ysofx   = p(5);
% Yeo     = k(13);
Yex     = p(6);


%FIXED
qm = 0.0067;
Yeo = 0.8163;
qEpmax = 0.635;
Kep = 0.5034;
Kis = 2.2729;
% Kec = 0.4267;
Yso = 0.4684;
Yse = 0.1955;
Ks = 0.045;
Kie = 1.6129;

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
dxdt(1) = X*(mu-D);                 % dXdt
dxdt(2) = D*(Sin-S)-qS*X;           % dSdt
dxdt(3) = qE*X-D*E;                 % dEdt
dxdt(4) = klao2*(osat-O)-qO*X-D*O;  % dOdt
dxdt(5) = F;


dxdt = dxdt';



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%