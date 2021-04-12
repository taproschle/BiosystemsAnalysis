%TODO: put the units and the description of each parameter and equation
function dxdt = anane_model(t, x)

X = x(1);
S = x(2);
A = x(3);
DOTa = x(4);
DOT = x(5);
V = x(6);

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
%Si     = 0; %g/L concentración de entrada

if t < 10
    F = 0;
    Si = 0;
else
    F = 0.035;
    Si = 300;
end

% Algebraic equations
qs = qS_max/(1 + A/Kia)*(S/(S + Ks)); 
qsof = Pa_max*qs/(qs + Kap);
qsox = (qs - qsof)*DOT/(DOT + Ko);
qsA = qa_max/(1 + qs/Kis)*A/(A + Ksa);
pa = qsof*Yas;
qa = pa - qsA; % Just for plot 

mu = (qsox - qm)*Yem + qsof*Yxsof + qsA*Yxa;
qo = (qsox - qm)*Yos + qsA*Yoa; % Para DOT

% Falta definir la ley de control (flujo de entrada)
% Tambien falta el KLa y H para el DOT 
% Falta el tau
% ODEs 
dxdt = zeros(6,1);
dxdt(1) = -F/V*X + mu*X;
dxdt(2) = F/V*(Si - S) - qs*X;
dxdt(3) = -F/V*A + qsA*X;
% Definicion parcial
DOT_ast = 99;
KLa = 220;
H = 10000;
Kp = 3600 / 35;
dxdt(4) = KLa*(DOT_ast - DOTa) - qo*X*H;
dxdt(5) = Kp*(DOTa - DOT);
dxdt(6) = F;
end





