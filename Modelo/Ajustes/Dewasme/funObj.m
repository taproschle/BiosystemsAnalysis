function f = funObj(k)

data = load('data.csv');

texp    = data(:,1)';
yexp    = data(:,2:5);

% Fixed parameters
muset   = 0.11;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.035;
Ko      = 0.0001;
v       = [muset X0 V0 Sin klao2 osat Ko];

% % Adjusted parameters (overflow)
% Kie     = 10.00;
% Yes     = 0.480;
% Kec     = 0.100;
% Ysofx   = 0.020;
% Yoe     = 1.104;
% Yxe     = 0.720;
% qOmax   = 0.256;
% Yosof   = 0.000;
% kof     = [Kie Yes Kec Ysofx Yoe Yxe qOmax Yosof];

% Initial conditions
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;
y0 = [X0 S0 E0 O0 V0];

tspan = texp;
fun = @(t,y) dewasme_unified(t,y,v,k);
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'NonNegative',[1,2,3,4,5]);

ran = randi(3);
if ran == 1
    points = ".";
elseif ran == 2
    points = "..";
else
    points = "...";
end
disp('Iterating'+" "+points)

[~ , Y] = ode15s(fun,tspan,y0,options);

wt = [1 10 1 1];

if ~ isequal(size(Y(:,1:4)),size(yexp))
    n = 1e10*ones(1,4);
    f = dot(wt,n);
else
    nX = norm(Y(:,1) - yexp(:,1));
    nS = norm(Y(:,2) - yexp(:,2));
    nE = norm(Y(:,3) - yexp(:,3));
    nO = norm(Y(:,4) - yexp(:,4));
    n = [nX nS nE nO].^2;
    f = dot(wt,n);
end
end