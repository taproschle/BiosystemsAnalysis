function f = funObj(k)

data = load('data.csv');

texp    = data(:,1)';
yexp    = data(:,2:4);

% Fixed parameters
muset   = 0.11;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.035;
Ko      = 0.0001;
v       = [muset X0 V0 Sin klao2 osat Ko];

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

wt = [1 1 1];

if ~ isequal(size(Y(:,1:3)),size(yexp))
    f = 1e20;
else
    nX = (Y(:,1) -yexp(:,1))./max(yexp(:,1));
    nS = (Y(:,2) - yexp(:,2))./max(yexp(:,2));
    nE = (Y(:,3) - yexp(:,3))./max(yexp(:,3));
    n = sum([nX nS nE].^2);
    f = dot(wt,n);
end
end