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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kL,k0,kU,opts_SSm,texp,ydata,x0,solver_ODE,opts_ODE,T,U] =...
    load_problem

%Lower, initial and upper values for each parameter in SSm:%

%  P = [ p1 p2    p3    p4    p5 p6 p7 p8 p9 p10 p11 p12 p13 ]
kL   = [0.1 0.1 0 0.1 0.1];
k0   = [3.9960 0.296 0.0051 1.7531 0.3497 ];
kU   = [5 1 0.1 5 1];

%SSm options:
opts_SSm.maxeval      = 1000;
opts_SSm.local.n1     = 100;
opts_SSm.maxtime      = 100;
opts_SSm.strategy     = 3;
opts_SSm.local.solver = 'fminsearch';
opts_SSm.local.finish = 'fminsearch';
opts_SSm.combination  = 1;  

%Experimental Data:
% [Time X S P]
data  = table2array(readtable('data.csv')); %edit here
data = data(2:length(data),1:5);
[~,n] = size(data); 
texp  = data(:,1)';
ydata = data(:,2:n);

%Initial conditions for integration:
%      X S E O 
x0 = [4.13 0.001 4.1 0.0024]; 

%ODE options:
solver_ODE = 'ode23';
opts_ODE   = odeset('RelTol',1e-4,'AbsTol',1e-7);

%Threshold for correlations (any couple of parameters with a correlation
%higher than this value will be fixed):
T = 0.95;

%Threshold for iterations (criterion III):
U = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%