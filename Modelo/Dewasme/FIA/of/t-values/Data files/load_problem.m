%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load_problem
% Define all the needed data and options to solve the problem
%
% OUTPUTS:
% kL            Row vector with the parameter lower bounds for SSm.
% k0            Row vector with the parameter initial values for SSm.
% kU            Row vector with the parameter upper bounds for SSm.
% opts_SSm      Options for the SSm procedure
% texp          Experimental times for the integration.
% ydata         Experimental data. Be consistent with obj_var.
% x0            Initial state variable vaues for the integration.
% solver_ODE    Solver for the dynamic model integration
% opts_ODE      Options for the dynamic model integration (odeset)
% T             Threshold for correlation values.
%
% Benjamín J. Sánchez: 2014-07-17
% Javiera Pérez: 17-05-2018
% Cristobal Torrealba: 30-06-2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kL,k0,kU,opts_SSm,texp,ydata,x0,solver_ODE,opts_ODE,int_time] = load_problem

%SSm options:
opts_SSm.maxeval      = 1000;
opts_SSm.local.n2     = 1;
opts_SSm.local.n1     = 2;
opts_SSm.maxtime      = 50;
opts_SSm.strategy     = 3;
opts_SSm.local.solver = 'fmincon';
opts_SSm.local.finish = 'fmincon';
opts_SSm.combination  = 1;  
opts_SSm.local.tol    = 1;  
opts_SSm.local.iterprint    = 1;  
opts_SSm.iterprint    = 1;  
%% Experimental Data:

%Lower, initial and upper values for each parameter in SSm:%

kB =  [ 1e-2, 1;     % Ks
        1e-2, 100;   % qSmax
        1e-2, 2;     % Ysoxx
        1e-5, 0.5;   % Yos
        1e-2, 10;    % Kie
        1e-2, 1;     % Yes
        1e-4, 1e-2;  % Kec
        1e-2, 0.5;   % Ysofx
        1e-2, 2;     % Yoe
        1e-2, 10;    % Yex
        1e-2, 2;     % qOmax
        1e-4, 1e-1]; % Yosof
    
k0      = [ 0.1000   % Ks
            3.5000   % qSmax
            0.4900   % Ysoxx
            0.3968   % Yos
            10.000   % Kie
            0.4800   % Yes
            0.1000   % Kec
            0.0200   % Ysofx
            1.1040   % Yoe
            0.7200   % Yxe
            0.2560   % qOmax
            100e-6]; % Yosof
        
kL   = kB(:,1);
kU   = kB(:,2);


%Temperature and integration time information
data = load('data.csv');
[~,n] = size(data); 
int_time  = data(:,1);
y_data = data(:,2:n);
%% Hippo stuff
data  = [int_time y_data]; 

[~,n] = size(data); 
texp  = data(:,1);
ydata = data(:,2:n);

%Initial conditions for integration:
X0 = 4.125;
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;
V0 = 0.3;
x0 = [X0 S0 E0 O0 V0];

%ODE options:
solver_ODE = 'ode45';
opts_ODE   = odeset('RelTol',1e-3,'AbsTol',1e-3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%