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
Kio     = 4;
v       = [muset X0 V0 Sin klao2 osat Ko Kio];

% % Adjusted parameters (overflow)
% Kie     = 5;
% Yes     = 0.667;
% Kec     = 0.05;
% qEmax   = 0.2;
% Ysofx   = 0.15;
% Yoe     = 1.067;
% Yxe     = 0.4;
% qOmax   = 0.43;
% Kio     = 4;
% kof     = [Kie Yes Kec qEmax Ysofx Yoe Yxe qOmax Kio];

% Initial conditions
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;
y0 = [X0 S0 E0 O0 V0];

% kID = ["Ks","qSmax","Ysoxx","qm","Yos","Kie","Yes","Kec",...
%         "qEmax","Ysofx","Yoe","Yxe","qOmax"];
% tab = table(kID',k');
% disp(tab)

tspan = texp;
fun = @(t,y) xu_unified(t,y,v,k);
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