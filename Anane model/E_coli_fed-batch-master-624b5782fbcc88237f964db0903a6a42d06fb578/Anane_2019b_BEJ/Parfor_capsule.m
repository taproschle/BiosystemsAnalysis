function [MC_results_raw,Fval] = Parfor_capsule(k1,data)

number_startspoints = 1;
%------------------- Fermentation conditions --------------------

Si              = 300;                            % concentration of glucose feed 300 g/L
mu_set          = 0.222;                          % set spec. growth rate during exp. feed
DOTstar         = 99;
Kla             = 220;
tau             = 25;
data.u          = [Si mu_set DOTstar Kla tau];
data.Fe0        = 0.145*60/1000;                  % L/h; initial value of exponential feed.
data.freq       = 1;
data.wt         = [1  1  1  0.1];                 % weighting factors for X, S, A, DoT
data.off_idx    = [1    19   20   38  39  57];    % indices of offline data X = [1 19], S = [20 38], A = [39 57]

%------------------------- Initialization --------------------------------

data.y0         = [2  0.17 4.94 0.0129 98 98 0 ];   %Initial conditions [V0, X0, S0,  A0  DOT0 DOTm F0]. 


                 %[Kap,1     Ksa,2   Ks,3    Kia,4   Kis,5  pAmax,6  qAmax,7   qm,8 qSmax,9  Yas,10  Yoa,11 Yxa,12  Yem,13  Yos,14  Yxsof,15  ] 
ParIni          = [0.4050   0.0730   0.1000  4.0000  8.000  0.1500  0.1700   0.0500  0.8000  0.9472  1.0670  0.8000  0.5100   1.0670 0.150];
% ParIni          = [1.4613   0.2360   0.0287  0.8340  0.200  0.0723  0.5146   0.0324  1.2251  0.9016  0.6947  0.4237  0.7000   1.3701 0.162];

lb              = [0.1000   0.0050   0.0050  0.1000  0.100  0.0500  0.0300   0.0100  0.400   0.2000  0.5000  0.2000  0.3500   0.2000  0.050];        
ub              = [0.9000   0.3000   0.5100  20.000  20.00  0.5000  0.6000   0.1000  1.500   0.9900  1.8000  0.9900  0.7000   2.0000  0.300];      
data.Par_ref    = [0.4050   0.0730   0.1000  4.0000  8.000  0.1500  0.1700   0.0500  0.8000  0.9472  1.0670  0.8000  0.5100   1.0670 0.150];
% nPar            = length(ParIni);
%------------ Normalize Parameters to the interval [0 1] -----------------

Par_norm    = ParIni./data.Par_ref;
lb_norm     = lb./data.Par_ref; 
ub_norm     = ub./data.Par_ref;

%------------- Parameter Estimation using fmincon ------------------------
    

    load([pwd '/Results/MonteCarlo_data.mat'],'MonteCarlo_data')
    data.Offlinedata_vec(:,2)  =  MonteCarlo_data(k1,1:data.off_idx(end))';
    data.Offlinedata_vec(:,1)  =  MonteCarlo_data(end,1:data.off_idx(end))'; %#ok<*COLND>

    data.RawDOT(:,2)  =  MonteCarlo_data(k1,(data.off_idx(end)+1):end)';
    data.RawDOT(:,1)  =  MonteCarlo_data(end,(data.off_idx(end)+1):end)';

    data.optroute = 0;
    
    %% use global optimization algorithm MS or GS to run a global search with each MC data.
    opts = optimoptions(@fmincon,'Algorithm','interior-point','OptimalityTolerance',1e-8,'UseParallel',false,'StepTolerance',1e-8);

    problem = createOptimProblem('fmincon','objective',@(Par_norm)Objective_function(Par_norm,data),'x0',Par_norm,'lb',lb_norm,'ub',ub_norm,'options',opts);

    MS = MultiStart;

    [Pmin,Fval] = run(MS,problem,number_startspoints);
   
    par_All = Pmin.*data.Par_ref;

    MC_results_raw = par_All;
    
    save([pwd '/Results/RawMC_results/RawMC_' num2str(k1) '_' num2str(Fval) '.mat'],'MC_results_raw')

end

