function PFR_STR_v8
% This code simulates a plug flow reactor (PFR) with a recycling stream from an STR.
% The PFR is modelled as a series of tiny CSTRs, each representing a single
% plug. The ode15s-solver is used to solve differential equations of a
% simple e. coli model for each plug/CSTR.
% Adaptation options:
% Total volume and geometry of the pipe, number of plugs, (initial
% conditions)


% Authors: Annina Sawatzki, Emmanuel Anane, M Nicolas Cruz B

% May 2018
%% for plotting
BW = 1.5;
FS = 20;

%% Initial values in STR

Inits(1) = 1.85;             % Biomass concentration, X [g/L]
Inits(2) = 1;                % Substrate concentration, S [g/L]
Inits(3) = 100;              % Dissolved oxygen, DOT [%]
Inits(4) = 0;                % Acetate concentration, A [g/L]
Inits(5) = 10;               % Volume [L]

data.num_var = length(Inits);            % number of variables

t_start = 0;                 % start of simulation [h]
t_end   = 50*90/3600;         % end of simulation [h]
t_step  = 1/3600;             % iteration step

%% Parameters:

%Parameter values from cultivations of E. coli published in Anane et. al 2017 (Biochem. Eng. J)

data.Par(01) = 0.349;       % Kap, monod-type saturation constant, intracellular acetate production (g/(g*h))
data.Par(02) = 0.020;       % Ksa, affinity constant, acetate consumption (g/L))
data.Par(03) = 0.019;      % Ks, affinity constant, substrate consumption (gglu./L)
data.Par(04) = 1.111;       % Kis, inhibition constant, inhib. of ace.uptake by glucose (g/L)
data.Par(05) = 287;         % Kla, volumetric mass transfer coefficient (per hour)
data.Par(06) = 0.203;       % pAmax, max spec acetate production rate (g_ace/(gx.h)))
data.Par(07) = 0.106;       % qAmax, max spec acetate consumption rate (g_ace/(gx.h)))
data.Par(08) = 0.013;      % qm, spec maintenance coefficient (g_glu/(gx.h))
data.Par(09) = 0.635;       % qSmax, max spec glucose uptake rate (g_glu./gx.h)
data.Par(10) = 0.827;      % Yas, yield of intracellular acetate on substrate (g.ace/g.glu) intracellular. How much of the total overflow flux is converted to acetate. The remaining goes to PPP and mixed acid pathways. 
data.Par(11) = 0.611;       % Yxa, yield of biomass on acetate
data.Par(12) = 0.546;       % Yem, yield exclusive maintenance (Yem = Yxs-qm)
data.Par(13) = 0.206;      % Yxsof, yield exclusive maintenance (Yem = Yxs-qm)
data.Par(14) = 0.250;       % Ypx, yield of product on biomass
data.Par(15) = 0.460;       % Yoa, yield, mass of oxygen consumed per gram of acetate oxidised
data.Par(16) = 0.479;       % Yos, yield, mass of oxygen consumed per gram of glucose oxidised
data.Par(17) = 0.0001;      % Ko

%% Variables
% PFR:
Vtotal = 1.8;                                          % total volume of the PFR [L]
N_pfr = 100;                                           % number of 'mini-STRs' in the PFR
V_pfr = Vtotal/N_pfr;                                  % volume per 'mini-STRs' [L]
r = (25*10^(-3))/2;                                    % [m]; diameter of pipe = 25 mm
A = pi * r^2;                                          % cross section of the pipe [m^2]
L_r = (V_pfr*10^(-3))/A;                               % length of each reactor [m]

data.N = N_pfr + 1;                                    % total number of reactors (PFR + STR)

% Point of feed addition: 
V_feedadd = 0.15;                                      % [L]
N_feedadd = round(V_feedadd/(Vtotal/N_pfr))+1;         % 'Mini-STR' which is closest to the calculated volume; 1 is added as for the STR (N=1), the PFR thus starts at N=2

% Feed
data.Xin = 0;                                          % biomass concentration in feed
data.Sin = 250;                                        % substrate concentration in feed
data.DOTin = 0;                                        % DOT in feed
data.Ain = 0;                                          % acetate concentration in feed
data.F = zeros(1,data.N);
data.F(N_feedadd) = 0.7;                               % Feed rate when feed is constant

% Flow rates, etc.
data.R  = 5 * ones(1,data.N);                           % flow rate in PFR between the STR and the feed inlet [L/h]
tau     = 3600 * V_pfr/data.R(N_feedadd+2);                % residence time for mini-STR at feed addition [s]

% Sampling points
V_where_sp = [0.24, 0.48, 0.72, 0.96, 1.2] + 0.15;     % Volume of the PFR where a sample is taken
V_sp = 0.00001;                                        % Sampling rate [L/h]
N_sp = round(V_where_sp./(Vtotal/N_pfr))+1;            % 'Mini-STR' which are closest to the calculated volume for 
data.Spl = zeros(1,data.N);
for k_sample = 1:length(N_sp)
    data.Spl(1,N_sp(k_sample)) = V_sp;
end

%% Initialize
% Vector of initial conditions for the mini-STRs
  Inits_pfr = [0, 0, 0, 0, V_pfr];
% Inits_pfr = [X  S DOT A  V]

y0 = [Inits, repmat(Inits_pfr,1,(data.N-1))];           % vector of initial conditions for all parts of the system

%% Numerical method

tspan = t_start:t_step:t_end;

% Preallocating the matrices
X = zeros(length(tspan),data.N);
S = zeros(length(tspan),data.N);
DOT = zeros(length(tspan),data.N);
A = zeros(length(tspan),data.N);
V = zeros(length(tspan),data.N);

X_str = zeros(length(tspan),1);
S_str = zeros(length(tspan),1);
DOT_str = zeros(length(tspan),1);
A_str = zeros(length(tspan),1);
V_str = zeros(length(tspan),1);

X_pfr = zeros(length(tspan),data.N-1);
S_pfr = zeros(length(tspan),data.N-1);
DOT_pfr = zeros(length(tspan),data.N-1);
A_pfr = zeros(length(tspan),data.N-1);
V_pfr = zeros(length(tspan),data.N-1);

options = odeset('NonNegative',1:length(y0));

[t,y] = ode15s(@ecoli,tspan,y0,options,data);    % the set of differential equations is solved for P_all(n_P)

X = y(:,1:data.num_var:end);       % chooses values of X
S = y(:,2:data.num_var:end);       % chooses values of S
DOT = y(:,3:data.num_var:end);     % chooses values of DOT
A = y(:,4:data.num_var:end);       % chooses values of E
V = y(:,5:data.num_var:end);       % chooses values of V

% STR
X_str = X(:,1);
S_str = S(:,1);
DOT_str = DOT(:,1);
A_str = A(:,1);
V_str = V(:,1);

% PFR
X_pfr = X(:,2:end);
S_pfr = S(:,2:end);
DOT_pfr = DOT(:,2:end);
A_pfr = A(:,2:end);
V_pfr = V(:,2:end);

%% Plotting

%% Plot of variable values over time in each reactor
%%{

figure
subplot(3,2,1)
plot(t*3600,X_str)
ylabel('X [g/L]')
xlabel('Time [s]')
title(['Concentrations in the STR at \tau = ', num2str(tau),' s'])
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

subplot(3,2,2)
plot(t*3600,S_str)
ylabel('S [g/L]')
xlabel('Time [s]')
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

subplot(3,2,3)
plot(t*3600,DOT_str)
ylabel('DOT [%]')
xlabel('Time [s]')
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

subplot(3,2,4)
plot(t*3600,A_str)
ylabel('A [g/L]')
xlabel('Time [s]')
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

subplot(3,2,5)
plot(t*3600,V_str)
ylabel('V [L]')
xlabel('Time [s]')
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 


%% 3D-plot of length (X-axis), time (Y-axis) and variable value (Z-axis) for the PFR
%%{

L_all = L_r*(1:data.N-1);
% view coordinates
az = -45;
el = 45;

% Plot of X

figure
M = t;
Y = L_all;
Z = X_pfr';
mesh(M,Y,Z)
xlabel('Time [h]')
ylabel('Length [m]')
zlabel('X_{PFR} [g/L]')
h = colorbar;
ylabel(h, 'Biomass concentration [g/L]')
view(az, el)
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 


% Plot of S
figure
M = t;
Y = L_all;
Z = S_pfr';
mesh(M,Y,Z)
xlabel('Time [h]')
ylabel('Length [m]')
zlabel('S_{PFR} [g/L]')
h = colorbar;
ylabel(h, 'Substrate concentration [g/L]')
view(az, el)
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

% Plot of DOT
figure
M = t;
Y = L_all;
Z = DOT_pfr';
mesh(M,Y,Z)
xlabel('Time [h]')
ylabel('Length [m]')
zlabel('DOT_{PFR} [g/L]')
h = colorbar;
ylabel(h, 'DOT [%]')
view(az, el)
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

% Plot of A
figure
M = t;
Y = L_all;
Z = A_pfr';
mesh(M,Y,Z)
xlabel('Time [h]')
ylabel('Length [m]')
zlabel('A_{PFR} [g/L]')
h = colorbar;
ylabel(h, 'Acetate concentration [g/L]')
view(az, el)
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

%}

end

function dydt = ecoli(t,y,data)

    num_var = data.num_var;
    N = data.N;
    R = data.R;                 % 'base' flow rate
    Spl = data.Spl;             % flow rate of sampling
    F = data.F;                 % flow rate of feed addition
    Xin = data.Xin;
    Sin = data.Sin;
    DOTin = data.DOTin;
    Ain = data.Ain;

    X = y(1:num_var:end);       % choses values of X
    S = y(2:num_var:end);       % choses values of S
    DOT = y(3:num_var:end);     % choses values of DOT
    A = y(4:num_var:end);       % choses values of A
    V = y(5:num_var:end);       % choses values of V_str
    
    % Parameters: 
    Kap     = data.Par(1);      % monod-type saturation constant, intracellular acetate production (g/(g.h))
    Ksa     = data.Par(2);      % affinity constant, acetate consumption (g/L)
    Ks      = data.Par(3);      % affinity constant, substrate consumption (gglu./L)
    Kis     = data.Par(4);      % inhibition constant, inhib. of ace.uptake by glucose(g/L)
    Kla     = data.Par(5);      % volumetric mass transfer coefficient (per hour)
    pAmax   = data.Par(6);      % max spec acetate production rate (g_ace/(gx.h)))
    qAmax   = data.Par(7);      % max spec acetate consumption rate (g_ace/(gx.h)))
    qm      = data.Par(8);      % spec maintenance coefficient (g_glu/(gx.h))
    qSmax   = data.Par(9);      % max spec glucose uptake rate (g_glu./gx.h)
    Yas     = data.Par(10);     % yield of intracellular acetate on substrate (g.ace/g.glu) intracellular. How much of the total overflow flux is converted to acetate. The remaining goes to PPP and mixed acid pathways. 
    Yxa     = data.Par(11);     % yield of biomass on acetate
    Yem     = data.Par(12);     % yield exclusive maintenance (Yem = Yxs-qm)
    Yxsof   = data.Par(13);     % biomass yield from the overflow route
    Ypx     = data.Par(14);     % yield of product on biomass
    Yoa     = data.Par(15);     % yield, mass of oxygen consumed per gram of acetate oxidised
    Yos     = data.Par(16);     % yield, mass of oxygen consumed per gram of glucose oxidised    
    Ko      = data.Par(17);
    
    % Constants
	Cs      = 0.391;
    Cx      = 0.488;
    H       = 14000;
    DOTstar = 90;                                   % based on values of empty wells 

    % Algebraic equations:
    qS      = (qSmax.*S)./(S+Ks);
    qSox    = (qS-(pAmax*(qS./(qS+Kap)))).*(DOT./(DOT+Ko));
    qSof    = qS-qSox;            
    pA      = qSof*Yas;      
    qSan    = (qSox-qm)*Yem*Cx/Cs;  
    qsA     = (qAmax./(1+(qS./Kis))).*(A./(A+Ksa));
    qA      = pA-qsA;
    mu      = (qSox-qm)*Yem + qsA*Yxa + (qSof-pA)*Yxsof;
    qO      = Yos*(qSox-qSan)+qsA*Yoa;
    
    %% Differential equations:
    
    % STR (j=1)
    R(1) = R(N)+F(1)-Spl(1);
    dydt(1:num_var,1) = [-X(1)*R(1)/V(1) + X(N)*R(N)/V(1) + Xin*F(1)/V(1) + mu(1)*X(1);
                         -S(1)*R(1)/V(1) + S(N)*R(N)/V(1) + Sin*F(1)/V(1) - qS(1)*X(1);
                         -DOT(1)*R(1)/V(1) + DOT(N)*R(N)/V(1) + DOTin*F(1)/V(1) + Kla(1)*(DOTstar-DOT(1)) - qO(1)*X(1)*H;
                         -A(1)*R(1)/V(1) + A(N)*R(N)/V(1) + Ain*F(1)/V(1) + qA(1)*X(1);
                         sum(F)-sum(Spl)];
    
    % PFR (j=2 to j=N)
    for j = 2:N      
        R(j) = R(j-1)+F(j)-Spl(j);
        dydt(((j*num_var)-(num_var-1)):(j*num_var),1) =  [-X(j)*R(j)/V(j) + X(j-1)*R(j-1)/V(j) + Xin*F(j)/V(j) + mu(j)*X(j);
                                                         -S(j)*R(j)/V(j) + S(j-1)*R(j-1)/V(j) + Sin*F(j)/V(j) - qS(j)*X(j);
                                                         -DOT(j)*R(j)/V(j) + DOT(j-1)*R(j-1)/V(j) + DOTin*F(j)/V(j) - qO(j)*X(j)*H;           % no aeration
                                                         -A(j)*R(j)/V(j) + A(j-1)*R(j-1)/V(j) + Ain*F(j)/V(j) + qA(j)*X(j);
                                                         0];
    end
    
end