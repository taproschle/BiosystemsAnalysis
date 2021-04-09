load mbp6.mat; load Pulse_time_2XP.mat

data.PE_status          = 'no';  %'yes' for PE, 'no' for simulation and plotting
data.Metabolic_states   = 'yes'; %Plot metabolic states in simulations?
%% Parameter Estimation conditions
Global_PE               = 'yes'; %'yes' for global PE, 'no' for local solution
data.freq               = 1;
data.wt                 = [1.5  1  1  0.2 1];    %weighting factors for X, S, A, DoT, P,
data.off_idx            = [1  27  28  54  55  81  82  105];         %indices of offline data X = [1 8], A = [9 16], P = [17 24]
number_startspoints     = 20;
MS_runs                 = 3;


%%  Experimental data: collect data of all triplicate bioreactors of this condition
data.Offlinedata    = [];
Rx_num              = [6 14 22];
Rxtor_num                = 14;

data.tspan               = Bioreactor(Rxtor_num).DOT(25:5040,1);
data.DOT                 = Bioreactor(Rxtor_num).DOT(25:5040,2);

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

%% Pick out event times
[~, FBinit]                = min(abs(data.tspan-6.2956));  %these time points are valid for all 24 bioreactors 

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

data.y0      = [2.90404 0   0.105  data.DOT(data.kla_change(6))  0];
% data.Konst  = [0.0603    0.07588   5.646  8.869  0.06635  0.7803  0.5773 466];
              %Ksa       Ks         Kip    Kis    qm     Yas     Ypx     Kla

% ParIni     = [0.860  0.0135  0.0048  0.878   1.580   0.3107   0.507     0.041   1.9489   0.246   0.221    0.449   0.1978     0.542    0.957   1.169   300  200  300  350 400]; % Original guess

             %[Kap,1   Ksa,2    Ks,3    Kip,4   Kis,5     pAmax,6  qAmax,7   qm,8   qSmax,9  Yas,10    Yxa,11   Yem,12  Yxsof,13   Ypx,14   Yoa,16  Yos,17  Kla1, 18--21 ] 
% ParIni     = [1.0587  0.03805  0.0598  1.5197   3.62811   0.7843   0.1701   0.0607  1.9370   0.5761    0.5304   0.5121  0.4134     0.3934   1.2718  1.276   643];%   573  701 432  346]; % PE1
% ParIni      = [0.1169  0.7074  0.2834  1.634   0.6893   0.4155  0.561    1.217   1.441]; %PE2
% lb          = [0.100   0.100   0.100   1.200   0.100    0.400   0.100    0.100   0.10];%  150  250  300 300];        
% ub          = [1.500   0.850   0.650   3.500   0.800    0.800   0.900    1.500   1.80];%  800  800  800 800];      
               %[Kap     pAmax   qAmax  qSmax   Yxa     Yem     Yxsof    Yoa     Yos] 
load Pmin.mat
ParIni = Pmin(1,:);

              %[Kap     Ksa     Ks    Kip    Kis     pAmax    qAmax  qm        qSmax   Yas     Yxa     Yem     Yxsof   Ypx    Yoa     Yos kla] 
% ParIni      = [0.1169  0.0603 0.0758 5.646  8.869  0.7074     0.2834   0.06635 1.634 0.9703  0.1893  0.5155    0.16136 0.5773 1.217 1.441 10  475  523]; %PE2
% ParIni     = [0.2580 0.0498 0.0563 4.3444 6.3751 0.5898 0.1033 0.0421 1.3066 0.5356 0.1503 0.5362 0.8843 0.3157 2.3628 2.4953 10  475  523]; 
% ParIni     = [0.497412727793123,0.0200617814679632,0.0398372570756547,3.95678736196967,1.30414757914609,0.328112624523344,0.117964671922563,0.0337145574372697,2.32914847523908,0.732498809511869,0.212408160528296,0.671372358963235,0.216756510439669,0.349664078949334,2.76349978014211,2.35002387046722,47.9376562590565,433.979797008276,530.537479823478];
% ParIni  =[0.525551375561463,0.0156761248999657,0.0430379787950411,4.36269273591996,6.24286366785587,0.376187287508553,0.111833763526454,0.0374949417929959,1.68538607890167,0.758704611700084,0.205913189677266,0.655320185544821,0.219244606525458,0.300953806707762,2.85846503654370,2.58259495008371,52.1864621116029,475.359208261873,603.498314299023];
% ParIni = [0.533382029756698,0.0158845029427112,0.0438420698679996,4.60962794504628,6.47543271580201,0.402806347489739,0.115290022926551,0.0371856577811353,1.72062594405969,0.768674607759207,0.209631188090950,0.656725736654459,0.220054165002372,0.300959872881055,2.90386344087627,2.63630903247540,52.2566792637010,474.599880150807,608.829338714042];
% ParIni = [0.595398473982403,0.0783896689127496,0.0501671460037222,4.55215595460179,6.93942272368842,0.431636787831069,0.119937312116475,0.0393919736901581,1.80564386753668,0.741104630981802,0.426836492660621,0.694135237854981,0.225214679389269,0.300078288709650,3.18549191005680,3.28346777930610,50.1414097646159,478.771266937239,668.820823738470];
% ParIni = [0.584524348767554,0.0930693993576722,0.0523657967050496,4.44390572391095,7.02285014731777,0.370601708140951,0.105413709125829,0.0393751255127791,1.56840412547423,0.739479089277735,0.460461787897354,0.685930898527602,0.237580935912097,0.287850585769475,3.19810796340619,3.32294876530519,49.9944260735085,486.518478516583,672.185287205395];

lb          = [0.100   0.010    0.01   1.00    1.00   0.100   0.08     0.01   0.800   0.100   0.100    0.45   0.10    0.100    0.800   0.80 30 300 400];%  150  250  300 300];        
ub          = [1.00   0.100    0.080  10.00   10.0   0.850   0.650     0.080   3.500  0.990   0.800    0.800  0.300   0.690    4.500   4.50  60 600 800];%  800  800  800 800];      


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
        data.optroute = 0;
        opts = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',1e-8,'UseParallel',false,'StepTolerance',1e-8,'Display','iter');
        [x2, fval] = fmincon(@(Par_norm)ObjFn_mbp6(Par_norm,data),Par_norm,[],[],[],[],lb_norm,ub_norm,[],opts);

        Pmin_final = x2.*data.Par_ref;

        save('Pmin_final.mat','Pmin_final')   
    end

end
