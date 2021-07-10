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
ModelName='Anane';

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

%Here, we define model parameters as th{i} for a total of i = 1..15 regression parameters.
for i=1:15
    syms(sprintf('th%d',i),'real');
end

%Additionally, we define variables t symbolic values representing time and model inputs.
syms t real

%Once symbolic parameters and model states are defined, we build arrays to vectorize parameters and model states.
statesSym=[x1 x2 x3 x4 x5]';
thetaSym=[th1 th2 th5 th6 th7 th8 th9 th10 th11 th13 th14 th15]';

%Define the amount of model states and parameters
nx=length(statesSym);
nth=length(thetaSym);

%Define initial conditions vector as a symbolic variable associated to numerical values.
thIC=sym(x0);
x0Sym=thIC';

%Define the amount of simulations performed in the Monte Carlo procedure
NExp=100; 

%% Nominal parameter values (In this example, we use those in Zenteno et al. (2010))

thetaNom = [ 1.9097     % qSmax     / th1
             0.242      % qEcmax    / th2
             0.045      % Ks        / th5
             0.4267     % Kec       / th6
             0.5034     % Kep       / th7
             1.6129     % Kie       / th8
             2.2729     % Kis       / th9
             0.2032     % YSoxX     / th10
             0.1949     % YSofX     / th11
             0.8163     % Yeo       / th13
             0.5485     % Yex       / th14
             0.1955];   % Yse       / th15
                        % X         / x1
                        % S         / x2
                        % E         / x3
                        % O         / x4
                        % V         / x5
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
th12    = 0.4684; % Yso     Reparam1
th3     = 0.635;  % qEpmax
th4     = 0.0067; % qm      Reparam2

%NOTE: This is the ODE system representing the Zenteno et al. (2010) model using no abreviatons for terms such as specific growth rate (mu).

qS      = th1*x2/(x2+th5)*1/(1+x3/th8);
qSof    = th3*qS/(qS+th7);
qSox    = (qS-qSof)*x4/(x4+Ko);
qEp     = qSof*th15;
qEc     = th2*x3/(x3+th6)*1/(1+qS/th9);
qO      = (qSox - th4)*th12 + qEc*th13;
mu      = (qSox-th4)*th10 + qSof*th11 + qEc*th14;
F0      = muset*(X0*V0)/(th10*Sin);
F       = F0*exp(muset*t);
D       = F/x5;

Xdot =  [ x1*(mu-D);
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
