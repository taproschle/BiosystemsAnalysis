%Otro intento
function dydt= intento_anane(t, y, Par,u)

%global flag1 flag2  flag3

%Emmanuel Anane, Diana C Lopez C, Peter Neubauer, M Nicolas Cruz Bournazou
%Chair of Bioprocess Enginering,
%Institute fur Biotechnologie, TU-Berlin.
%28/01/2017

%This function is the set of ODEs for a fed-batch culture of e coli. The script that 
%calls it implements a fed-batch process of wild type E. coli fermentation
%without recombinant protein production.
%%
V       = y(1);           % volume
X       = y(2);           % biomass
S       = y(3);           % substrate
A       = y(4);           % acetate
DOT     = y(5);           % dissolved oxygen, no probe response
DOTm    = y(6);           % dissolved oxygen, with probe response time
%F       = y(7);          % feed
%%  Fermentation constants
Cs = 0.391;%              %carbon content of glucose (ref:enfors)
Cx = 0.488;               %carbon content of biomas(e coli)  
H  = 14000;               %Henry's law constant to convert dissolved oxygen from mol to % DO
%% INPUTS
Si          = u(1);                     
mufeed      = u(2);
DOTstar     = u(3);     
Kla         = u(4);
tau         = u(5);
%% Parameters

Kap     = Par(1);   %monod-type saturation constant, intracellular acetate production
Ksa     = Par(2);   %affinity constant, acetate consumption (0.05 g/L)
Ko      = Par(3);   %affinity constant, oxygen consumption (molO./L)
Ks      = Par(4);   %affinity constant, substrate consumption (0.05 gglu./L)
Kia     = Par(5);   %inhibition constant, inhib. of glucose uptake by acetate
Kis     = Par(6);   %inhibition constant, inhib. of ace.uptake by glucose(4g/L)
pAmax   = Par(7);   %max spec acetate production rate (0.15 g_ace/(gx.h)))
qAmax   = Par(8);   %max spec acetate consumption rate (0.15 g_ace/(gx.h)))
qm      = Par(9);   %spec maintenance coefficient (0.04 g_glu/(gx.h))
qSmax   = Par(10);  %max spec glucose uptake rate (1.3 g_glu./gx.h)
Yas     = Par(11);  %yield of acetate on substrate (g.ace/g.glu)
Yoa     = Par(12);  %yield of oxygen on acetate g/g)
Yxa     = Par(13);  %yield of biomass on acetate
Yem     = Par(14);  %yield exclusive maintenance (Yem = Yxs-qm)
Yos     = Par(15);  %yield of oxygen on glucose (g/g)
Yxsof   = Par(16);  %biomass yield from the overflow route

%% Explicit algebraic equations (See text in paper for explanations)

qS      = (qSmax./(1+(A/Kia)))*(S/(S+Ks)); % Check
qSof    = pAmax*(qS/(qS+Kap));             % Check
pA      = qSof*Yas;                        % Check
qSox    = (qS-qSof)*(DOT/(DOT+Ko));        % Check      
qSan    = (qSox-qm)*Yem*Cx/Cs;             % 
qsA     = (qAmax./(1+(qS/Kis)))*(A/(A+Ksa));
qA      = pA-qsA;                            
mu      = (qSox-qm)*Yem + qsA*Yxa + qSof*Yxsof;
qO      = Yos*(qSox-qSan)+qsA*Yoa;
Kp      = (1/tau)*3600;

% Acá va el tema de los regimenes de alimentacion 
if t < 11.44                    %batch phase
    F = 0;                      %no feeding
elseif t < 16.3                 %exponential fed-batch phase, non-production.
    mu_set = 0.3;
    F = mu_set/(Yxa*Si)*2*0.17*exp(mu_set*t);
else           %constant feed fed-batch phase, protein induction
    F = 0.0128;
end
%% ODEs
dydt = zeros(6,1);
dydt(1,1) = F;                          %dV/dt
dydt(2,1) = X*(mu - (F/V));             %dX/dt
dydt(3,1) = F/V*(Si - S) - (qS*X);      %dS/dt
dydt(4,1) = qA*X-(F/V)*A;               %dA/dt
dydt(5,1) = Kla*(DOTstar-DOT)-qO*X*H;   %dDOT/dt
dydt(6,1) = Kp*(DOT-DOTm);              %dDOTa/dt

end