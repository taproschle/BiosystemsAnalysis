%{
Framework for Parameter estimation, ill-conditioning diagnosis,
regularization and output uncertainty analaysis of dynamic growth models,
as presented in the paper ('   ').
%}
%% Step 1: Parameter Estimation (PE) using a multistart approach 

load RawDOT; load Offlinedata_vec; load Offlinedata.mat

mkdir('Results')

data.RawDOT           = RawDOT;
data.Offlinedata_vec  = Offlinedata_vec;
data.Offlinedata      = Offlinedata;
data.freq             = 1;
data.wt               = [1  1  1  0.1];    %weighting factors for X, S, A, DoT
data.off_idx          = [1    19   20   38  39  57]; %indices of offline data X = [1 19], S = [20 38], A = [39 57]

MS.Random_seeding     = 1;            % 1 for 'YES', 0 for 'NO'
MS.startpoints        = 10;           % number of random seeds
MS.runs               = 1;

data.Tikh_Reg         = 0;             % 1 for 'Tikhonov regularization', 0 for 'No regularization'

              %[Kap,1     Ksa,2   Ks,3    Kia,4   Kis,5  pAmax,6  qAmax,7   qm,8 qSmax,9  Yas,10  Yoa,11 Yxa,12  Yem,13   Yos,14   Yxsof,15] 
Para.Ini     = [0.4050   0.0730   0.1000  4.0000  8.000  0.1500  0.1700   0.0500  0.8000  0.9472  1.0670  0.8000  0.5100   1.0670  0.150];
Para.lb      = [0.1000   0.0050   0.0050  0.1000  0.1000 0.0500  0.0300   0.0100  0.400   0.4000  0.5000  0.2000  0.3500   0.2000  0.050];        
Para.ub      = [2.0000   0.3000   0.5100  5.0000  10.000 0.5000  0.6000   0.1000  1.500   0.9900  1.8000  0.9900  0.7000   2.0000  0.300];      


[Pmin_all, Fval_all] = PE_multistart(data,MS,Para);

%-----------Plot profiles with optimal parameter values-------------------
% Pmin = [0.5088	0.0128	0.0381	1.2602	1.8383	0.2286	0.1148	0.0133	0.635	0.8938	0.5221	0.5794	0.5321	1.5722	0.229]; %best Pmin
Pmin = Pmin_all(1,:);
S_e_coli(Pmin,data)
nPar = length(Pmin);

%% Step 2: Calculate Absolute Sensitivity Matrix at Optimal Parameter Values
Rank        = 1;
repeat      = 0;
tspan       = data.RawDOT(:,1);
count       = 1;

while Rank ~= length(Pmin)
    
    Sens_calc      = 'yes';
    Par            = Pmin;
    y0             = [2 0.17 4.94 0.0129 98 98 0];  
    y0jac          = zeros(1,nPar*length(y0));

    y0             = horzcat(y0,y0jac);

    [~,~,Sens]     = S_e_coli_output(tspan,Par,y0,Sens_calc);

    save([pwd '/Results/Sens ' datestr(now) '.mat'],'Sens')

 %% Step 3  Singular Value Decomposition & Rank Calculation -------------

    CondNMaxCrit        = 1000; % Condition Number Grah's threshold (Lopez et al 2016)
    CollinIndexMaxCrit  = 15;   % Collinearity Index Brun's threshold between 10-15 (Brun 2002)

    Rank   = Rank_calculation(Sens,CondNMaxCrit,CollinIndexMaxCrit);

    %% Step 4: Uncertainty analysis for well-posed problems [Rank = length(Par)]
    if Rank == length(Pmin) && repeat == 0

       %----calculate covariance matrix, Std dev. and confidence intervals----
        dof         = length(Sens)-nPar;            % degree of feedom
        s_squar     = Fval(1)./dof;                 % estimator variance
        Cov_matrix  = s_squar*inv(Sens'*Sens) ;     % covariance matrix of parameters
        Sigma_Par   = sqrt(diag(Cov_matrix))';      % standard deviation of parameter estimates
        alfa        = 0.025;                        % significance level,
        tcr         = tinv((1-alfa),dof);           % critical t-dist value at alfa
        lb_est      = Pmin-Sigma_Par*tcr;           % -95% confidence intervals
        ub_est      = Pmin+Sigma_Par*tcr;           % +95% confidence intervals

        %----calculate correlation matrix.-------------------

        Corr_mat = zeros(nPar,nPar);
        for i6 = 1:nPar
            for i7 = 1:nPar
                Corr_mat(i6,i7) = Cov_matrix(i6,i7)/sqrt(Sigma_Par(i6)^2*Sigma_Par(i7)^2);
            end
        end

       %-----generate correlated samples from the interval [lb ub]------
       N_samples     = 500;  
       Samples       = lhsdesign(N_samples,nPar); 
       New_samples   =  bsxfun(@times,Samples,ub_est-lb_est);
       Ind_samples   =  bsxfun(@plus,New_samples,lb_est);
       Corr_samples  =  imancon(Corr_mat,Ind_samples);                        

       %----run model simulations using correlated samples, plot and save results --------------
       par_label = {'K_{ap}','K_{sa}','K_s','K_{ia}', 'K_{is}', 'p_{Amax}', 'q_{Amax}', 'q_m', 'q_{Smax}', 'Y_{as}', 'Y_{oa}', 'Y_{xa}', 'Y_{em}', 'Y_{os}', 'Y_{xsof}'};
       MC_raw  = 'no';
       Sens_calc = 'no';
       TableOutput_error = Output_uncertainty_analysis(N_samples,par_label,tspan,Ind_samples,Corr_samples,MC_raw,Sens_calc);

    else            
    %% Step 5: Regularization and Monte Carlo Analysis for ill-conditioned Problems [Rank < length(Par)]
    
       if Rank == length(Pmin)
           break
       end
       
       Tune_vector  = [1e2 1e4  1e6 10  1]; % vector for tuning the regularization
       repeat       = 1;
       
%-----------Generate MC data by Bootstrap sampling-------------------%
%{       
       N_samples     = 500;  

       MonteCarlo_data             = Boostrap_sampling(Pmin,data,N_samples);

        
       save([pwd '/Results/MonteCarlo_data.mat'],'MonteCarlo_data')
%}     
%-----Generate MC data by random sampling using Std dev and mean----------%
       
        N_samples = 20;                                        % number of Monte Carlo Samples to generate.
        DOT_std_dev                 = RawDOT(:,2).*0.05;        % since dissolved oxygen is a high precision measurement, its std dev is assumed to be 5 of the measured values%
        Standard_devs               = vertcat(Offlinedata_vec(:,3),DOT_std_dev);
        Measurements                = vertcat(Offlinedata_vec(:,2),RawDOT(:,2));

        MonteCarlo_data             = MC_Sampling(N_samples,Measurements,Standard_devs);
        MonteCarlo_data(end+1,:)    = vertcat(Offlinedata_vec(:,1),RawDOT(:,1)); % for the time
        save([pwd '/Results/MonteCarlo_data.mat'],'MonteCarlo_data')

%----Run Regularized Parameter Estimation with each Monte Carlo dataset---%
        mkdir('Results','RawMC_results')
        data.Tikh_Reg      = 1;             % 1 for 'Tikhonov regularization', 0 for 'No regularization'
        data.Tikh_lambda   = Tune_vector(count);
        nPar               = 15;
        MC_results_raw     = zeros(nPar,N_samples);
        Fval               = zeros(1,N_samples);

        parfor k1          = 1: N_samples;

               [MC_results_raw(:,k1),Fval(:,k1)] = Parfor_capsule(k1,data);

        end

        [C, I]          = sort(Fval,'ascend');

        MC_results_sorted    = MC_results_raw(:,I);    %sort according Fval

        par_All             = MC_results_sorted;
        save([pwd '/Results/MC_results_All ' datestr(now) '.mat'],'MC_results_sorted')


        %---------  Plot results of MC parameter estimation in 3D----------%

        par_label = {'K_{ap}','K_{sa}','K_s','K_{ia}', 'K_{is}', 'p_{Amax}', 'q_{Amax}', 'q_m', 'q_{Smax}', 'Y_{as}', 'Y_{oa}', 'Y_{xa}', 'Y_{em}', 'Y_{os}', 'Y_{xsof}'};

        alpha     = 0.05;       % confidence level
        Ranking   = 2;          % no ranking

        PE_results_3Dplot(par_All,par_label,alpha,Ranking)
        
        %------- Calculate mean and standard deviation of MC_estimates----%
        avg_Par = mean(par_All,2);
        Pmin    = avg_Par;   
        save([pwd '/Results/MC_results_avg ' datestr(now) '.mat'],'avg_Par')

    end
 count = count + 1;  
end

%% Step 5.1: Uncertainty Analysis of Model Outputs with MC_data
if (repeat)
    MC_raw         = 'yes';
    Sens_calc      = 'no';
    Corr_samples   = par_All;
    par_label = {'Kap','Ksa','Ks','Kia', 'Kis', 'pAmax', 'qAmax', 'qm', 'qSmax', 'Yas', 'Yoa', 'Yxa', 'Yem', 'Yos', 'Yxsof'};


    TableOutput_error = Output_uncertainty_analysis(N_samples,par_label,tspan,[],Corr_samples,MC_raw,Sens_calc);
end


%% Framework completed