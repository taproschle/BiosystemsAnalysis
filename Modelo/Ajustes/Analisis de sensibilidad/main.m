
clc,clear
%% ANANE
load data.csv
texp    = data(:,1)';
yexp    = data(:,2:5);

tsim    = texp(end);

load 'kAn.mat'

X0      = 4.125;
V0      = 0.3;
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;
y0 = [X0 S0 E0 O0 V0];
tspan   = [0 tsim];
fun = @anane_unified;
flag = 3;
NumEqs = [1 2 3 4];
figure(1)
ksensibilidad(NumEqs,k,y0,tspan,fun,flag) 
figure(2)
identifica(NumEqs,k,y0,tspan,fun,0.95)


%% XU

clc,clear

load data.csv
texp    = data(:,1)';
yexp    = data(:,2:5);
load 'kXu.mat'

tsim    = texp(end);


X0      = 4.125;
V0      = 0.3;
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;
y0 = [X0 S0 E0 O0 V0];
tspan   = [0 tsim];
fun = @xu_unified;
flag = 3;
NumEqs = [1 2 3 4];
figure(1)
ksensibilidad(NumEqs,k,y0,tspan,fun,flag) 
figure(2)
identifica(NumEqs,k,y0,tspan,fun,0.95)


%% DEWASME


clc,clear

load data.csv
texp    = data(:,1)';
yexp    = data(:,2:5);
load 'kDew.mat'

tsim    = texp(end);


X0      = 4.125;
V0      = 0.3;
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;
y0 = [X0 S0 E0 O0 V0];
tspan   = [0 tsim];
fun = @dewasme_unified;
flag = 3;
NumEqs = [1 2 3 4];
figure(1)
ksensibilidad(NumEqs,k,y0,tspan,fun,flag) 
figure(2)
identifica(NumEqs,k,y0,tspan,fun,0.95)
