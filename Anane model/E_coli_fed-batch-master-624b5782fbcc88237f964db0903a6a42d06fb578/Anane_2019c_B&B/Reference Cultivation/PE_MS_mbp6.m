load mbp6.mat;

data.PE_status          = 'no';  %'yes' for PE, 'no' for simulation and plotting
data.Metabolic_states   = 'no'; %Plot metabolic states in simulations?
%% Specifications for PE and simulations
Global_PE               = 'yes'; %'yes' for global PE, 'no' for local solution
% data.freq               = 1;
data.wt                 = [1  1  1  0.2 1];    %weighting factors for X, S, A, DoT, P,
data.off_idx            = [1  30  31  60  61  90  91  111];  %indices of offline data X = [1 9], S = [10 23]...
number_startspoints     = 20;
MS_runs                 = 5;
Rx_num                  = [2 10 18];   %reference cultivation
%% Experimental data: collect data of all triplicate bioreactors of this condition

Rxtor_num               = 10;
tspanOrig              = Bioreactor(Rxtor_num).DOT(25:4970,1);
DOTOrig                = Bioreactor(Rxtor_num).DOT(25:4970,2);

data.Offlinedata = [];

for k0 = 1:length(Rx_num)
    
Plot_data(k0).biomass =  Bioreactor(Rx_num(k0)).OD600(2:end,:);

data.Offlinedata         = [data.Offlinedata;...
                           Plot_data(k0).biomass];
                    
end

for k0 = 1:length(Rx_num)
Plot_data(k0).Glu =  Bioreactor(Rx_num(k0)).Glucose(2:end,:);

data.Offlinedata        = [data.Offlinedata;...
                            Plot_data(k0).Glu];
end
                            
for k0 = 1:length(Rx_num)
Plot_data(k0).Ace =  Bioreactor(Rx_num(k0)).Acetate(2:end,:);

data.Offlinedata        = [data.Offlinedata;...
                            Plot_data(k0).Ace];
end

for k0 = 1:length(Rx_num)
Plot_data(k0).Pro =  Bioreactor(Rx_num(k0)).Product(2:end,:);
    
data.Offlinedata        = [data.Offlinedata;...
                            Plot_data(k0).Pro];
end                            

data.Plot_data = Plot_data;
%% Extract event times during cultivation                            
   
[~, FBinit]                = min(abs(tspanOrig-5.35122));  %these time points are valid for all 24 bioreactors 
data.tspan                 = tspanOrig(FBinit:end); % re-define simulation time to start at start of fed-batch
data.DOT                   = DOTOrig(FBinit:end);


data.kla_change(1)         = 1;

[~, data.kla_change(2)]    = min(abs(data.tspan-8.3036));

[~, TInduce]               = min(abs(data.tspan-8.93333));  %these time points are valid for all 24 bioreactors 


data.kla_change(3)         = TInduce; %put this in the kla change vector, but it is not a kla change event, it is induction time

%% for fed-batch phase only

%% Inputs (Experimental conditions)

Si                  = 200;
F                   = 0;
V                   = 0.01;
mu_set              = 0.35;
pH                  = 6.87;

data.u              = [Si  V F mu_set pH];

data.U0             = 0.2;
%% Initialization

% data.y0          = [0.037 5   0  data.DOT(1)  0  0];

               %[Ksa     Ks      Kip      Kis     qm      Yas  Klas 
data.Konst = [0.0134  0.00487  1.531  1.58026  0.02594  0.9856];
%ParIni     = [0.860    0.0135  0.0048  0.878 1.5802   0.3107  0.157   0.04094  1.9489 0.24566  0.22186  0.44906 0.19788 0.54244  0.95752 1.16299  1.30   7.5  300  200  300  350 400];
%ParIni = [0.764004379380186,0.0104118355258620,0.00975992700636704,0.937559450012777,2.77654817628846,0.402634179087940,0.143283204501285,0.0398242959045009,1.95568913063593,0.655024407664758,0.447793108151367,0.461023247812516,0.229133925298310,0.475987561448710,0.937541754758169,1.15465570968839,1.85142584542458,6.67350119531736,504.696600966012,275.294832819845,312.830104467162,342.071724316039,434.542535786551];
% ParIni  = [0.745286122146964,0.00601641975468218,0.0302942320306567,2.76806113384588,4.53574010181202,0.608909749347876,0.260326924729023,0.0286461283677788,1.80902550648700,0.811766572406549,0.458089415061957,0.518713341671113,0.188844530864908,0.402684122412182,0.979875350839119,1.18218219624048,2.00718175660914,7.83697114769062,538.196062131934,192.556774783998,271.857591426644,322.822899651135,654.963859962606];
%ParIni     = [1.021    0.0472   0.0334 2.9640 7.875   0.644   0.329   0.0202   2.0289 0.9017   0.4079   0.4627   0.3586  0.4694  1.2761   1.1433  1.1823 7.161 474 342  437  463 651];
%ParIni     = [1.021   0.644   0.329   2.0289    0.4079    0.4627   0.3586  0.4694  1.2761   1.1433  1.1823 7.161];
% ParIni     = [0.735953118517400,0.615927843922152,0.329795153683726,2.08430470379284,0.661943697496441,0.422566917299501,0.430882092629384,0.424285612611608,1.01737938814766,1.32707658179553,1.12726870623932,6.82196122263109 474 342  437  463 651 780];
% ParIni     = [0.6428 0.6453  0.3299  2.0468  0.7524  0.4342  0.4304  0.3918 0.9977 1.3801 1.0385 6.4228 515 294 369 482 589 780];
% ParIni     = [1.05561059766142,0.584191027013514,0.314835585139339,2.00737020473958,0.735743266229447,0.446672273177751,0.406217795212519,0.424148449093979,1.10697681254694,1.43169319807844,1.02202593939202,6.14675831747499 564 780];
          %[Kap  Ksa     Ks      Kip      Kis  pAmax   qAmax   qm   qSmax   Yas Yxa     Yem    Yxsof    Ypx     Yoa     Yos     U1  U2  kla1 kla2 kla3 kla4 kla5] 

ParIni = [1.1173 0.0134  0.00487  1.531  1.58026  0.5603  0.2599  0.02594  1.9401 0.9856 0.7231  0.4489  0.3988  0.5189  0.9445  1.5003  1.0202  5.5059  601.36  601.4108];  

% load Pmin.mat
% ParIni = Pmin(1,:);
% ParIni = [0.650348418130514,0.529116860120405,0.376729039158638,2.01639520824455,0.802244613148980,0.435175491809657,0.317913050221099,0.403121973042734,1.00284291045356,1.36808607228498,1.03203457942530,6.05555420628879,479.640695348824,273.058765814636,396.805061616981,455.356221008856,580.819446758352,745.131522542162];
lb         = [0.500  0.001   0.001 1.00 1.00 0.100   0.100    0.002  1.20  0.574  0.100    0.4      0.10    0.200   0.500    0.500   0.10   1.00  200  250];%  200 250  300 400];        
ub         = [1.850  0.080   0.010 20.0 20.0 0.800   0.800    0.080  2.80   0.99  0.999    0.7      0.500   0.700   1.500     1.800  5.00   10.0  800   800];%  800 800 800 1000];
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
    if strcmp(Global_PE,'yes')     

        Par_min    = zeros(MS_runs,length(ParIni));                      %Preallocations
        Fval       = zeros(MS_runs,1); 

        for k3 = 1:MS_runs
            data.optroute = 0;
            opts = optimoptions(@fmincon,'Algorithm','interior-point','OptimalityTolerance',1e-8,'UseParallel',true,'StepTolerance',1e-8);

            problem = createOptimProblem('fmincon','objective',@(Par_norm)ObjFn_mbp6(Par_norm,data),'x0',Par_norm,'lb',lb_norm,'ub',ub_norm,'options',opts);

    %         MS = MultiStart;
            GS = GlobalSearch;
    %         [x,fval] = run(MS,problem,number_startspoints);
            [x,fval] = run(GS,problem);

            Par_min(k3,:)    = x.*data.Par_ref;
            Fval(k3)         = fval; 

        end
        [Fval_ordered, I]          = sort(Fval,'ascend'); % Rank the sum of squared errors (SSE) or residuals in ascending order,

        Pmin            = Par_min(I,:);      % select the optimal parameter values that gave the minimum residual.

        save('Pmin.mat','Pmin')
        save('Fval_ordered.mat','Fval_ordered')

    else
        %% Local PE with lsqnonlin
        data.optroute = 0;
        opts = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',1e-8,'UseParallel',false,'StepTolerance',1e-8,'Display','iter');
        [x2, fval] = fmincon(@(Par_norm)ObjFn_mbp6(Par_norm,data),Par_norm,[],[],[],[],lb_norm,ub_norm,[],opts);

        Pmin_final = x2.*data.Par_ref;

        save('Pmin_final.mat','Pmin_final')   
    end
end

