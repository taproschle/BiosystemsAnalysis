function dydt= fn_e_coli(t, dydt, Par,u)

global flag1 flag2  flag3
%Emmanuel Anane, Nicolas Cruz B. M.
%28/07/2016

%this function is the set of ODEs for a fed-batch culture. dydt(1) = dV/dt,
%dydt(2) = dx/dt, dydt(3) = dS/dt, dydt(4) = dP/dt, dydt(5) = dF/dt. The script that 
%calls it implements a fed-batch process of e. coli fermentation, with
%induction at a later stage in the fermentation.
%%
V = dydt(1);           %volume
X = dydt(2);           %biomass
S = dydt(3);           %substrate
A = dydt(4);           %acetate
DOT = dydt(5);         %dissolved oxdydtgen, no probe response
DOTm = dydt(6);        %dissolved oxdydtgen, with probe response time
F = dydt(7);           % feed
P = dydt(8);           %Intracellular rProtein product
%%  Fermentation constants
Cs = 0.391;%                    %carbon content of glucose (ref:enfors)
Cx = 0.488;                     %carbon content of biomas(e coli)  
H = 14000;                      %Henry's law constant to convert from mol to % DO
%% INPUTS
Si = u(1);
mufeed = u(2);
DOTstar = u(3);     
I = u(4);           %inducer concentration
Xind = u(5);        %biomass at time of induction
%%
% Parameters
Kap = Par(1);
Kaq = Par(2);
Ko = Par(3);
Ks = Par(4);
Ksq = Par(5);
mumax = Par(6);
pAmax = Par(7);
qAmax = Par(8);
qm = Par(9);
qSmax = Par(10);
Yaof = Par(11);
Yaresp = Par(12);
Yem = Par(13);
Yosresp = Par(14);
Kla = Par(15);
tau = Par(16);
k1   =  Par(16);
k2 = Par(17);
%%
%Explicit algebraic equations
mu = mumax*S/(S+Ks);
qS = qSmax*S/(S+Ks)*DOT/(DOT+Ko);
qSof = pAmax*qS/(qS+Kap)/Yaof;
pA = pAmax*qS/(qS+Kap);                 %acetate production
qSox = qS-qSof;                        
qSan =(qSox-qm)*Yem*Cx/Cs;      
qsA = qAmax*A/(A+Kaq).*(Ksq/(Ksq+qS));  %acetate consumption
qA = pA-qsA;                            %net acetate accumulation
qO = Yosresp*(qSox-qSan)+qsA*Yaresp;
Kp = (1/tau)*3600;
qp = (I/Xind)*(k1*mu + k2); 
%%
dydt = zeros(8,1);
dydt(1,1) = F;                          %dV/dt
dydt(2,1) = X*(mu - (F/V));             %dX/dt
dydt(3,1) = F/V*(Si - S) - (qS*X);      %dS/dt
dydt(4,1) = qA*X-(F/V)*A;               %dA/dt
dydt(5,1) = Kla*(DOTstar-DOT)-qO*X*H;   %dDOT/dt
dydt(6,1) = Kp*(DOT-DOTm);              %dissolved oxygen with probe response 
if flag1 == true                        %batch phase
    dydt(7,1) = 0;                      %no feeding
    dydt(8,1) = 0;                      %no product formation 
elseif flag2 == true   %exponential fed-batch phase, non-production.
    dydt(7,1) = mufeed*F;               %dF/dt
    dydt(8,1) = 0;                      %no product formation
elseif flag3 == true;  %constant feed fed-batch phase, protein induction
    dydt(7,1) = 0;   
    dydt(8,1) = -mu*P+qp;               %dP/dt, intracellular product
end

end


