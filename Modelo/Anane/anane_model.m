function dxdt = anane_model(t, x, u)

X = x(1);
S = x(2);
A = x(3);
DOTa = x(4);
DOT = x(5);
V = x(6);

X0 = u(1);
V0 = u(2);
Sin = u(3);
muset = u(4);


% Parameters
Kap    = 0.5052;
Ksa    = 0.0134;
Ko     = 0.0001;
Ks     = 0.0370;
Kia    = 1.2399;
Kis    = 2.1231;
Pa_max = 0.2268;
qa_max = 0.1148;
qm     = 0.0129;
qS_max = 0.6356;
Yas    = 0.9097;
Yoa    = 0.5440;
Yxa    = 0.5718;
Yem    = 0.5333;
Yos    = 1.5620;
Yxsof  = 0.2268;
% Parameters for DOTa ode
DOT_ast = 99;
KLa = 220;
H = 10000;
Kp = 3600 / 35;

% Algebraic equations
qs = qS_max/(1 + A/Kia)*(S/(S + Ks)); 
qsof = Pa_max*qs/(qs + Kap);
qsox = (qs - qsof)*DOT/(DOT + Ko);
qsA = qa_max/(1 + qs/Kis)*A/(A + Ksa);
pa = qsof*Yas;
qa = pa - qsA; % Just for plot 

mu = (qsox - qm)*Yem + qsof*Yxsof + qsA*Yxa;
qo = (qsox - qm)*Yos + qsA*Yoa; % Para DOT

%fprintf("El qsa es %f \n", qsA)

F = muset/(Yxa*Sin)*X0*V0*exp(muset*t);

%fprintf("Feed %f en t = %f \n", F, t)

if S <= 0; S = 0; end
% ODEs
qsA = qsA / 10;
dxdt = zeros(6,1);
dxdt(1) = -F/V*X + mu*X;
dxdt(2) = F/V*(Si - S) - qs*X;
dxdt(3) = -F/V*A + qsA*X;        % Acetato
dxdt(4) = KLa*(DOT_ast - DOTa) - qo*X*H;
dxdt(5) = Kp*(DOTa - DOT);
dxdt(6) = F;
end





