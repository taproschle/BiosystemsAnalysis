
function dx = model(t,x,k)

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

% State variables
X = x(1);   % Biomass (g/L)
S = x(2);   % Substrate (g/L)
E = x(3);   % Ethanol (g/L)
O = x(4);   % Dissolved Oxygen (g/L)
V = x(5);   % Volume (L)

% Fixed parameters
muset   = 0.11;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.035;
Ko      = 0.0001;

% Adjusted parameters (all cases)
Ks      = p(1);
qSmax   = p(2);
Ysoxx   = p(3);
Yso     = p(4);
Kio     = p(5);
Yse     = p(6);
Kec     = p(7);
Ysofx   = p(8);
Yeo     = p(9);
Yex     = p(10);
qOmax   = p(11);
Yosof   = p(12);

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
dx(1) = X*(mu-D);                 % dXdt
dx(2) = D*(Sin-S)-(qSox-qSof)*X;  % dSdt
dx(3) = (Yse*qSof-qE)*X-D*E;      % dEdt
dx(4) = klao2*(osat-O)-D*O-...    % dOdt
    (Yso*qSox+Yosof*qSof+Yeo*qE)*X;
dx(5) = F;                        % dVdt

dx = dx';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%