global  kx1 kx2 kx3 ks1 ks2 ke2 ke3 ko1 ko2 ko3 kc1 kc2 kc3 kos koe rs_max muo Ko  Ks  Ke  Kie  CO2s alfa beta gamma delta Henry Si P 
Si  = 500;                           % Feed glucose concentration, [g/L], R1
%                                         ==============  Open Loop
mu_set = 0.13;               % specific grow rate [1/h]
P = 1;                                  % Pressure [atm]

%  Yield coefficients original
% Biomass
kx1 = 0.49;                                         %  [g of X/g of S],   R1
kx2 = 0.05;                                         %  [g of X/g of S],   R2
kx3 = 0.72;                                         %  [g of X/g of E],   R3
%kx3 = 0.05;                                      % Carcamo [g of X/g of E],   R3

% Substrate
ks1 = 1;                                           %  [g of S/g of X],               R1
ks2 = 1;                                           %  [g of S/g of X],               R2

% Product
%ke2 = 0.48;                                %  [g of E/g of S],   R2
ke3 = 1;                                          %  [-],               R3
ke2 = 0.32;                                   % Carcamo   [g of E/g of S]

% Oxygen [Primer rendimiento super importante - Determina el overflow]
ko1 = 0.427;                                      %  [g of O2/g of S],  R1 
ko2 = 0;                                              %  [g of O2/g of S],  R2
ko3 = 1.104;                                       %  [g of O2/g of E],  R3
%ko3 = 0.03944;                              %  Carcamo [g of O2/g of E],  R3

% Carbon dioxide
kc1 = 0.5897;                              %  [g of CO2/g of S], R1
kc2 = 0.4621;                              %  [g of CO2/g of S], R2 
kc3 = 0.6249;                              %  [g of CO2/g of E], R3

kos = ko1;                                    %  [g of O2/g of S],  R1 
koe = ko3;                                   %  [g of O2/g of E],  R3

%  Kinetic parameters
rs_max = 3.5/4;                          % [g of S / g of X / h],  R1
muo = 0.256/2;                          % [g of O2/g of X/h], R1
mu_crit = muo*kx1/ko1;          % without oxygen limitations

% Equilibrium parameters
Ko  = 0.0001;                          % [g of O2/L],        R1
Ks  = 0.1;                                 % [g of S/L],         R1
Ke  = 0.1;                                 % [g of E/L],         R1
Kie  = 10;                                 % [g of E/L],         R1

% Henry Constants CO2 - O2
HenryCO2 = 0.7498;          %   [(L*atm)/gCO2] 
yCO2 = 0.5;
CO2s = yCO2 /HenryCO2; % 0.66684 [gCO2/L]
%CO2s =1.286;                    % [gCO2/L] Dewasme

% O2
Henry =26.409;                   %  [(L*atm)/gO2] 

% Control law
alfa = 0.034;                              % [1/h], R2
beta = 1.33;                                % [1/h], R2
gamma = 0.603;                       % [1/h], R2
delta = 0.89;                              % [-],   R3