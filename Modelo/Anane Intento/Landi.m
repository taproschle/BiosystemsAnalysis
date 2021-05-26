function [dydt] = Landi(t, x) %Poner paramtros
% A 0.16 grow rate

t0 = 10;

mu_max = 0.2;
Ks = 1.8;
beta = 0.89;
Y_XES = 0.68;
Y_XG1 = 0.57/0.69;
Y_XG2 = 0.14/0.92;
Y_EXP = 0.022/0.052;
kd = 0;
Gr = 450;
% ODEs

X = x(1); %Biomass
G = x(2); %Glucose
E = x(3); %Ethanol
V = x(4); %Volume

alpha = 0.5;
mu_set = alpha;

X0 = 3.63;
V0 = 2;
if t > 10
    mu_G = 0.18;
    Y_XG = Y_XG1;
    F0 = mu_set*X0*V0/(Y_XG*Gr);
else
    mu_G = 0.16*exp(-0.024*(t-t0));
    Y_XG = Y_XG2;
    F0 = mu_set*X0*V0/(Y_XG*Gr);
end


F = F0*exp(alpha*t);

dydt = ones(4,1);
%Biomass 
dydt(1) = -F/V*X + (mu_G + mu_max*E/(Ks + E)*(1-beta) - kd)*X;
%Glucose
dydt(2) = F/V*(Gr - G) - (mu_G/Y_XG)*X;
%Ethanol
dydt(3) = -F/V*E + (-mu_max*E/(Ks + E)/Y_XES*(1-beta) + mu_G*Y_XES)*X;
%Volume
dydt(4) = F;
%dydt(5) = (O - Oast) - kla*X;
end