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
% Kie     = 1.2399;
% pEmax   = 0.2268;
% Kep     = 0.5052;
% Yes     = 0.9097;
% Kec     = 0.0134;
% qEmax   = 0.1148;
% Kis     = 2.1231;
% Ysofx   = 0.2268;
% Yoe     = 0.5440;
% Yxe     = 0.5718;
% kof     = [Kie pEmax Kep Yes Kec qEmax Kis Ysofx Yoe Yxe];

% Initial conditions
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;
y0 = [X0 S0 E0 O0 V0];

tspan = texp;
fun = @(t,y) anane_unified(t,y,v,k);
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

wt = [1 1 1 10];

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