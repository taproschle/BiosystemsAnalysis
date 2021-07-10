
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
V = x(5);   % Volume (L)


% Fixed parameters
muset   = 0.197;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.035;
Ko      = 0.0001;

% Adjusted parameters
qSmax   = p(1);
Ysoxx   = p(2);
Yse     = p(3);
Ysofx   = p(4);



% FIXED

Yeo     = 1.998;
Yosof   = 0.0001;
Yso     = 0.3438;
Kio     = 4.7323;
Ks      = 0.1414;
Kec     = 0.0051;
Yex     = 1.7531;
qOmax   = 0.3497;


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
dxdt(1) = X*(mu-D);                 % dXdt
dxdt(2) = D*(Sin-S)-(qSox-qSof)*X;  % dSdt
dxdt(3) = (Yse*qSof-qE)*X-D*E;      % dEdt
dxdt(4) = klao2*(osat-O)-D*O - ...
    (Yso*qSox+Yosof*qSof+Yeo*qE)*X; % dOdt
dxdt(5) = F;                        % dVdt
                 

dxdt = dxdt';



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%