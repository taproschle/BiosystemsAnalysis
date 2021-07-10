%Adapted from https://sourceforge.net/projects/minimal-output-sets/
%by Cristóbal Torrealba 2021

%Application of the fast structural identifiability asessment applyed
%to the Zenteno et al. (2010) wine fermentation model. 
%Instructions for this example:

    % Run step 1, wait until all variables are generated. If wanted,
    % adjust the number of simulations (NExp); takes few seconds with
    % Zenteno´s model
    %
    % Run step 2, wait until all simulations are performed. If wanted,
    % adjust the bounds associated to random nominal parameter values
    % perturbation; takes few minutes with Zenteno´s model.
    %
    % Run step 3, wait until all graphics are generated. If wanted,
    % change measured states (sensors); takes a minute or so. 
    % Enjoy.

%% Initial procedures
clear ; clc ; close all
%Prepare computation for Zenteno model (NOTE: Parameter vector does not include initial conditions as unknown parameters!)
tic
ModelName='Xu';

%Initial conditions for integration
X0 = 4.125;
S0 = 0.001;
E0 = 4.104;
O0 = 0.007;
V0 = 0.3;
x0 = [X0 S0 E0 O0 V0];

%Here operational data corresponding to this example is loaded
% load('Operational_Data.mat')

%Define final integration time (Tf, in hours for this example) and
%number of integration steps (Nt) 
Tf=40;
Nt=100;

%% Symbolic Toolbox Variable defintions

%Here, we define model states as x{i} for a total of i = 1..5 model states.
for i=1:5
    syms(sprintf('x%d',i),'real');
end

%Here, we define model parameters as th{i} for a total of i = 1..12 regression parameters.
for i=1:14
    syms(sprintf('th%d',i),'real');
end

%Additionally, we define variables t symbolic values representing time and model inputs.
syms t real

%Once symbolic parameters and model states are defined, we build arrays to vectorize parameters and model states.
statesSym=[x1 x2 x3 x4 x5]';
thetaSym=[th1 th2 th3 th4 th5 th6 th7 th8 th9 th10 th11 th12 th14]';

%Define the amount of model states and parameters
nx=length(statesSym);
nth=length(thetaSym);

%Define initial conditions vector as a symbolic variable associated to numerical values.
thIC=sym(x0);
x0Sym=thIC';

%Define the amount of simulations performed in the Monte Carlo procedure
NExp=100; 

%% Nominal parameter values (In this example, we use those in Zenteno et al. (2010))

thetaNom = [ 1.6041   % qSmax  / th1
             0.1224   % qEcmax / th2
             0.2569   % qOmax  / th3
             0.0575   % qm     / th4
             0.1034   % Ks     / th5
             0.5019   % Kec    / th6
             4        % Kio    / th7
             7.4378   % Kie    / th8
             0.2014   % YSoxX  / th9
             0.6751   % YSofX  / th10
             0.8965   % Yso    / th11
             0.7607    %Yeo    / th12
             0.3652]; % Yse    / th14
                      % X      / x1
                      % S      / x2
                      % E      / x3
                      % O      / x4
                      % V      / x5
%% MODEL DEFINITION
%  Fixed parameters
muset   = 0.11;
X0      = 4.125;
V0      = 0.3;
Sin     = 450;
klao2   = 180*100;
osat    = 0.035;
Ko      = 0.0001;
a       = 100;
th13    = 0.155; % Yex Reparam1

%NOTE: This is the ODE system representing the Zenteno et al. (2010) model using no abreviatons for terms such as specific growth rate (mu).


qS      = th1*x2/(x2+th5)*1/(1+x3/th8);
qOs      = th3*x4/(x4+Ko)*1/(1+x3/th7);
qScrit  = ((qOs/th11)-th4*th9)/(1-th9);
% Usamos Smooth Max/Min (Max a+ / Min a-)
% smax = (x1*exp(a*x1) + x2*exp(a*x2))/(exp(a*x1) + exp(a*x2));
qSox    = (qS*exp(-a*qS)+qScrit*exp(-a*qScrit))/(exp(-a*qS)+exp(-a*qScrit));
qSof    = (0*exp(a*0)+(qS-qSox)*exp(a*(qS-qSox)))/(exp(a*0)+exp(a*(qS-qSox)));
qEp = (qSof-qSof*th10)*th14;
qeci    = th2*x3/(x3+th6);
qecj = (th3-qOs)*th12/(1-th13);
qEc      = (qeci*exp(-a*qeci)+qecj*exp(-a*qecj))/(exp(-a*qeci)+exp(-a*qecj));
qO = qOs + (qEc-qEc*th13)*th12;
mu      = (qSox-th4)*th9 + qSof*th10 + qEc*th13;
F0      = muset*(X0*V0)/(th9*Sin);
F       = F0*exp(muset*t);
D       = F/x5;

Xdot =  [   x1*(mu-D);
            D*(Sin-x2)-qS*x1;
            (qEp-qEc)*x1-D*x3;
            klao2*(osat-x4)-D*x4-qO*x1 - D*x4;
            F];

%% END MODEL DEFINITION AND AUXILIAR FUNCTION GENERATION
%(May take a few minutes to process)

%Here, gradients are calculated symbolically and automatically stored as separated Matlab scripts for later usage.
dfdxSym     = jacobian(Xdot,statesSym);
dfdthSym    = jacobian(Xdot,thetaSym);
dx0dthSym   = jacobian(x0Sym,thetaSym);

f = matlabFunction(Xdot,'vars',{t,statesSym,thetaSym},'file',ModelName);
dfdx = matlabFunction(dfdxSym,'vars',{t,statesSym,thetaSym},'file',['dfdx',ModelName]);
dfdth = matlabFunction(dfdthSym,'vars',{t,statesSym,thetaSym},'file',['dfdth',ModelName]);
IC = matlabFunction(x0Sym,'vars',{thetaSym},'file',['x0',ModelName]);
dICdth = matlabFunction(dx0dthSym,'vars',{statesSym,thetaSym},'file',['dICdth',ModelName]);
