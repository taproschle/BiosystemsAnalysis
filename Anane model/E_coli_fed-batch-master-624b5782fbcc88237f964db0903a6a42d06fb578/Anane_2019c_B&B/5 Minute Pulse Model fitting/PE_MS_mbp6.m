load mbp6.mat; load Pulse_time_12

data.PE_status          = 'no';  %'yes' for parameter estimation, 'no' for simulation and plotting
data.Metabolic_states   = 'yes'; %Plot metabolic states in simulations?
%% Parameter Estimation conditions
Global_PE               = 'yes'; %'yes' for global PE, 'no' for local solution
data.freq               = 1;
data.wt                 = [1  1  1  0.5 1];    %weighting factors for X, S, A, DoT and P respectively
data.off_idx            = [1  27  28  54  55  81  82  105];  %indices of offline data X = [1 27], S = [28 54], etc
number_startspoints     = 20;       % only used for multi-start global parameter estimation
MS_runs                 = 3;        %number of global PEs to run
%% Experimental data: collect data of all triplicate bioreactors of this condition

Rx_num = [3 11 19];         %minibioreactor numbers
data.tspan               = Bioreactor(Rx_num(2)).DOT(25:4970,1);
data.DOT                 = Bioreactor(Rx_num(2)).DOT(25:4970,2);


data.Offlinedata = [];

for k0 = 1:length(Rx_num)
    
Plot_data(k0).biomass =  Bioreactor(Rx_num(k0)).OD600(3:end,:);

data.Offlinedata         = [data.Offlinedata;...
                           Plot_data(k0).biomass];
                    
end

for k0 = 1:length(Rx_num)
Plot_data(k0).Glu =  Bioreactor(Rx_num(k0)).Glucose(3:end,:);

data.Offlinedata        = [data.Offlinedata;...
                            Plot_data(k0).Glu];
end
                            
for k0 = 1:length(Rx_num)
Plot_data(k0).Ace =  Bioreactor(Rx_num(k0)).Acetate(3:end,:);

data.Offlinedata        = [data.Offlinedata;...
                            Plot_data(k0).Ace];
end

for k0 = 1:length(Rx_num)
Plot_data(k0).Pro =  Bioreactor(Rx_num(k0)).Product(1:end,:);
    
data.Offlinedata        = [data.Offlinedata;...
                            Plot_data(k0).Pro];
end                            

data.Plot_data = Plot_data;
%% Extract event times during cultivation                            
[~, FBinit]                = min(abs(data.tspan-6.2122));  %these time points are valid for all 24 bioreactors 

data.kla_change(1)              = 1;   
[~, data.kla_change(2)]         = min(abs(data.tspan-2.957));   
[~, data.kla_change(3)]         = min(abs(data.tspan-3.819));   
[~, data.kla_change(4)]         = min(abs(data.tspan-4.448));   
[~, data.kla_change(5)]         = min(abs(data.tspan-4.943));   
data.kla_change(6)              = FBinit;   
data.Pulse_time                 = Pulse_time;     
%% Inputs (Experimental conditions)

Si                  = 200;
F                   = 0;
V                   = 0.01;
mu_set              = 0.35;

data.u              = [Si  V F mu_set];
   
%% Initialization

data.y0      = [2.854 0   0.105  data.DOT(data.kla_change(6))  0];



             %[Kap,1   Ksa,2    Ks,3    Kip,4   Kis,5     pAmax,6  qAmax,7   qm,8   qSmax,9  Yas,10    Yxa,11   Yem,12  Yxsof,13   Ypx,14   Yoa,16  Yos,17  Kla1, 18--21 ] 

% ParIni     = [0.860  0.0135  0.0048  0.878   1.580   0.3107   0.507     0.041   1.9489   0.246   0.221    0.449   0.1978     0.542    0.957   1.169   300  200  300  350 400]; % Original guess
% ParIni     = [1.0587  0.03805  0.0598  1.5197   3.62811   0.7843   0.1701   0.0607  1.9370   0.5761    0.5304   0.5121  0.4134     0.3934   1.2718  1.276   643];%   573  701 432  346]; % PE1
% ParIni     = [0.1169  0.7074  0.2834  1.634   0.6893   0.4155  0.561    1.217   1.441]; %PE2
% ParIni     = [0.1169  0.0603 0.0758 5.646  8.869  0.7074     0.2834   0.06635 1.634 0.9703  0.1893  0.5155    0.16136 0.5773 1.217 1.441]; %PE2
% ParIni     = [0.2580   0.0498  0.0563 4.3444 6.3751  0.5898   0.1033 0.0421   1.3066   0.5356  0.1503  0.5362  0.8843  0.3157 2.3628 2.4953 10  475  523]; 
% ParIni     = [0.438561881933890,0.0360483014893555,0.0693964395055391,3.90212591187590,7.44430183485035,0.592318634859392,0.142874905119685,0.0306353408203354,1.35209381741609,0.719932155241629,0.183249827419920,0.555027500499577,0.989766028960218,0.275418814482617,2.61131609204783,2.72406835319355,11.2684172163504,465.380056459624,685.331136070177];
% ParIni     = [0.451311091342338,0.0282932041348478,0.0487914269690589,3.54445301009638,6.20279839856812,0.609763409766010,0.112102431778099,0.0321909088465503,1.46212442881618,0.684798468406410,0.128503905116646,0.534248467060868,0.989880034980784,0.328730025972966,2.76490626635479,2.56928967788609,32.9282015142983,446.081749349159,240, 500, 654.548038187606];


load Pmin.mat           %load results of final PE
ParIni = Pmin(2,:);


lb          = [0.100   0.010    0.01   1.00    1.00   0.100  0.050     0.01   0.800   0.100   0.100    0.45   0.10    0.200    0.800   0.80 3  300  100 200 400];%  150  250  300 300];        
ub          = [1.00   0.100    0.080  10.00   10.0   0.850   0.650     0.080   3.500  0.990   0.800    0.800  0.990   0.690    3.500   3.50  60 600  400 600 800];%  800  800  800 800];      


data.Par_ref = ParIni;
%% Normalize Parameters to the interval [0 1]
Par_norm    = ParIni./data.Par_ref;
lb_norm     = lb./data.Par_ref; 
ub_norm     = ub./data.Par_ref;

%% Simulation XOR PE?
if strcmp(data.PE_status,'no')
    ObjFn_mbp6(Par_norm,data)
else

%% Global PE with fmincon and multistart (MS)/globalsearch (GS)
    if strcmp(Global_PE ,'yes')     
        
        Par_min    = zeros(MS_runs,length(ParIni));                      %Preallocations
        Fval       = zeros(MS_runs,1); 

        for k3 = 1:MS_runs
            data.optroute = 0;
            opts = optimoptions(@fmincon,'Algorithm','interior-point','OptimalityTolerance',1e-8,'UseParallel',true,'StepTolerance',1e-8);

            problem = createOptimProblem('fmincon','objective',@(Par_norm)ObjFn_mbp6(Par_norm,data),'x0',Par_norm,'lb',lb_norm,'ub',ub_norm,'options',opts);

%             MS = MultiStart;
            GS = GlobalSearch;

%             [x,fval] = run(MS,problem,number_startspoints);

            [x,fval] = run(GS,problem);

            Par_min(k3,:)    = x.*data.Par_ref;
            Fval(k3)         = fval; 

        end
        [Fval_sorted, I]     = sort(Fval,'ascend'); % Rank the sum of squared errors (SSE) or residuals in ascending order,

        Pmin                 = Par_min(I,:);      % select the optimal parameter values that gave the minimum residual.

        save('Pmin.mat','Pmin')
        save('Fval_sorted.mat','Fval_sorted')


    else
        %% Local PE with lsqnonlin
        data.optroute = 1;
        opts = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',1e-8,'UseParallel',false,'StepTolerance',1e-8,'Display','iter');
        [x2, fval] = fmincon(@(Par_norm)ObjFn_mbp6(Par_norm,data),Par_norm,[],[],[],[],lb_norm,ub_norm,[],opts);

        Pmin_final = x2.*data.Par_ref;

        save('Pmin_local.mat','Pmin_final')   
    end
end


