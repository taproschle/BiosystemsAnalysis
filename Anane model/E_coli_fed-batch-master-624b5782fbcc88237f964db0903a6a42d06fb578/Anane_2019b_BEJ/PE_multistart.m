function [Pmin, Fval] = PE_multistart(data,MS,Para)

Random_seeding   = MS.Random_seeding;
%% Experimental Inputs

Si              = 300;                        % concentration of glucose feed 300 g/L
mu_set          = 0.222;                        % set spec. growth rate during exp. feed
DOTstar         = 99;
Kla             = 220;
tau             = 25;
data.u          = [Si mu_set DOTstar Kla tau];

data.Fe0        = 0.145*60/1000;                  % L/h; initial value of exponential feed.


%% Initialization

data.y0    = [2  0.17 4.94 0.0129 98 98 0 ];   %Initial conditions
            %[V0, X0, S0,  A0  DOT0 DOTm F0]..   %Initial conditions

ParIni     = Para.Ini;

lb         = Para.lb;
ub         = Para.ub;

data.Par_ref = ParIni;
%% Normalize Parameters to the interval [0 1]
Par_norm    = ParIni./data.Par_ref;
lb_norm     = lb./data.Par_ref; 
ub_norm     = ub./data.Par_ref;
                  
if (Random_seeding)
    %% Global PE with fmincon and multistart
    startspoints     = MS.startpoints;
    MS_runs          = MS.runs;

    Par_min    = zeros(MS_runs,length(ParIni));                      %Preallocations
    Fval       = zeros(MS_runs,1);
    
    for k3 = 1:MS_runs
        data.optroute = 0;
        %     opts = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',1e-8,'UseParallel',true,'StepTolerance',1e-8);
        opts = optimoptions(@fmincon,'Algorithm','interior-point','OptimalityTolerance',1e-8,'UseParallel',true,'StepTolerance',1e-8);
        
        problem = createOptimProblem('fmincon','objective',@(Par_norm)Objective_function(Par_norm,data),'x0',Par_norm,'lb',lb_norm,'ub',ub_norm,'options',opts);
        
        MS = MultiStart;
        
        [x,fval] = run(MS,problem,startspoints);
        Par_min(k3,:)    = x.*data.Par_ref;
        Fval(k3)         = fval;
    end
%--- Rank the sum of squared errors (SSE) or residuals in ascending order--
    
    [Fval_sorted, I]    = sort(Fval,'ascend');
    
    Pmin                = Par_min(I,:);      %select the optimal parameter values that gave the minimum residual.
    save([pwd '/Results/Fval_sorted' datestr(now) '.mat'],'Fval_sorted')

else
    %% Local PE with fmincon
    data.optroute = 0;
    opts = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',1e-8,'UseParallel',true,'StepTolerance',1e-8,'Display','iter');
    [x, fval] = fmincon(@(Par_norm)Objective_function(Par_norm,data),Par_norm,[],[],[],[],lb_norm,ub_norm,[],opts);
    Pmin = x.*data.Par_ref;
    save([pwd '/Results/fval' datestr(now) '.mat'],'fval')

end
  save([pwd '/Results/Pmin' datestr(now) '.mat'],'Pmin')

end