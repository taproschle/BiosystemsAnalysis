function TableOutput_error = Output_uncertainty_analysis(N_samples,par_label,tspan,Ind_samples,Corr_samples,MC_raw,Sens_calc)

%{
Emmanuel Anane
Inst. of Biotechnology, Bioverfahrensteknik (BVT), TU-Berlin
25-Oct 2017
Finalized: 23-Jan 2018

This code computes random samples using the latin hypercube sampling
technique, taking into account correlation among the model parameters. The
random samples are drawn from the the input parameter space, defined by the
lower and upper bounds of parameters. The correlation matrix 
used in the script is the Rank correlation matrix, calculated from first
order linear approximation of the inverse of the Fisher information matrix.
MC-optimization.

The uncertainty analysis comprises using randomly sampled Par_values to
simulate the model. The variability of the output of the model in response
to perturbations in input parameter values is depicted qualitatively by 
plots showing one standard deviation from the mean and quantitatively by
the relative standard deviation of model outputs.

%}
%% Step 2.0: Preliminaries

% Par         = MC_results_summary.Mean_Par;              
% sigma       = MC_results_summary.Sigma_Par;
% Corr_mat    = MC_results_summary.Corr_mat;


FS  = 14;
PW  = 1.2;       %plotwidth: thickness of plotted line
BW  = 2;         %thickness of boundaries
col = [0.5 0.5 0.5 0.9];
lw  = 1.0;
%% Step 4.4: Plot Samples

if strcmp(MC_raw,'no')
    
        %----plot independent samples---

    CorS_fig = figure;
    lp = par_label;
    [h,ax,bax,P] = plotmatrix(Ind_samples) ;
    set(ax,'FontSize',8,'FontWeight','bold')

    for k1=1:length(par_label)
        ylabel(ax(k1,1),lp(k1))
        xlabel(ax(length(par_label),k1),lp(k1))
    end

    title(['Non-Correlated Latin Hypercube Samples, ',num2str(N_samples),' Samples'])
    set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

%     saveas(CorS_fig,[pwd '/Results/Corr_samples_plot' datestr(now) '.tif']);
%     save([pwd '/Results/Non_Corr_samples ' datestr(now) '.mat'],'Ind_samples')

    %-----------------------plot correlated samples--------------------
    CorS_fig = figure;
    lp = par_label;
    [h,ax,bax,P] = plotmatrix(Corr_samples) ;
    set(ax,'FontSize',8,'FontWeight','bold')

    for k1=1:length(par_label)
        ylabel(ax(k1,1),lp(k1))
        xlabel(ax(length(par_label),k1),lp(k1))
    end

    title(['Latin Hypercube Sampling with Correlation, ',num2str(N_samples),' Samples'])
    set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

%     saveas(CorS_fig,[pwd '/Results/Corr_samples_plot' datestr(now) '.tif']);
%     save([pwd '/Results/Corr_samples ' datestr(now) '.mat'],'Corr_samples')


else
    [m, N_samples] = size(par_All);
    CorS_fig2 = figure;
    lp = par_label;
    [h,ax,bax,P] = plotmatrix(Corr_samples) ;
    set(ax,'FontSize',8,'FontWeight','bold')

    for k1=1:m
        ylabel(ax(k1,1),lp(k1))
        xlabel(ax(m,k1),lp(k1))
    end

    title('Raw MonteCarlo PEs')
    set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

%     saveas(CorS_fig2,[pwd '/Results/RawMC_results' datestr(now) '.tif']);

end
%% Step 2.2: Run simulations using the sampled parameter values

X_output = zeros(length(tspan),N_samples); S_output = zeros(length(tspan),N_samples);
A_output = zeros(length(tspan),N_samples); DOT_output = zeros(length(tspan),N_samples);
T_par = zeros(length(tspan),N_samples); %only for the 'parallel pool'
y0 = [2 0.17 4.94 0.0129 98 98 0];   %Initial conditions
   %[V0, X0, S0, A0 DOT0 DOTm F0]


load T_dummy; load Y_dummy  %sometimes the parameter combination can't successfuly run a simulation
parfor j2               = 1:N_samples        % using parallel computing to improve computation speed
    
   try
    Par                 = Corr_samples(j2,:);
    [T_out, Y_output,~] = S_e_coli_output(tspan,Par,y0,Sens_calc);
    
    X_output(:,j2)      = Y_output(:,2);
    S_output(:,j2)      = Y_output(:,3);
    A_output(:,j2)      = Y_output(:,4);
    DOT_output(:,j2)    = Y_output(:,6);
    T_par(:,j2)         = T_out;    %not necessary for a normal 'for' loop
    
   catch
        T = T_dummy;                    % a few times the parameter combinations will let the solver crash, in this case use the optimal solution (at mean par values) as model outputs.
        Y = Y_dummy;                    % this occurs about 6-8 times out of 1000 samples, so the results are not affected.
        disp('Solver failed, using dummy outputs')
        X_output(:,j2)      = Y(:,2);
        S_output(:,j2)      = Y(:,3);
        A_output(:,j2)      = Y(:,4);
        DOT_output(:,j2)    = Y(:,6);
        T_par(:,j2)         = T;    %not necessary for a normal 'for' loop
   end   
end

T                       = T_par(:,1);  %not necessary for a normal 'for' loop


%% Step 2.3: Plot model outputs, calculate and plot standard deviation (or 10th and 90th percentiles)
Output1 = figure;
subplot(2,2,1)

X_mean = (mean(X_output,2))';
X_std  = (std(X_output,0,2))';
mseb(T',X_mean,X_std);
% hold on
% plot(Offlinedata(:,1),Offlinedata(:,2),'r*');hold off;
ylim([0 40]);xlim([2 35])
xlabel('Time(h)');ylabel('X (g/L)'); title('Biomass');
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

subplot(2,2,2)
S_mean = (mean(S_output,2))';
S_std  = (std(S_output,0,2))';
mseb(T',S_mean,S_std);
% hold on
% plot(Offlinedata(:,1),Offlinedata(:,4),'r*');hold off;
ylim([0 5]);xlim([2 35])
xlabel('Time(h)');ylabel('S (g/L)'); title('Substrate');
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

subplot(2,2,3)
A_mean = (mean(A_output,2))';
A_std  = (std(A_output,0,2))';
mseb(T',A_mean,A_std);
% hold on
% plot(Offlinedata(:,1),Offlinedata(:,4),'r*');hold off;
ylim([0 1]);xlim([2 35])
xlabel('Time(h)');ylabel('A (g/L)'); title('Acetate');
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

subplot(2,2,4)
DOT_mean = (mean(DOT_output,2))';
DOT_std  = (std(DOT_output,0,2))';
mseb(T',DOT_mean,DOT_std);
% hold on
% plot(RawDOT(:,1),RawDOT(:,7),'-.r')
ylim([0 100]);xlim([2 35])
xlabel('Time(h)');ylabel('DOT (g/L)'); title('DOT');
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

saveas(Output1,[pwd '/Results/Output_uncertainty_plot' datestr(now) '.tif']);
saveas(Output1,[pwd '/Results/Output_uncertainty_plot' datestr(now) '.fig']);

%% Calculate variations in model outputs (Prediction Error)

dev_X = (X_std./X_mean).*100;
[max_devX,I_X]  = max(dev_X);

dev_S = (S_std./S_mean).*100;
[max_devS,I_S]  =   max(dev_S);

dev_A = (A_std./A_mean).*100;
[max_devA,I_A]  = max(dev_A);

dev_DOT = (DOT_std./DOT_mean).*100;
[max_devDOT,I_DOT]  =  max(dev_DOT);

Output_error = num2cell([max_devX max_devS max_devA max_devDOT]);
Occurs_at    = num2cell([T(I_X) T(I_S) T(I_A)  T(I_DOT)]);

Output_label = {'Biomass','Glucose','Acetate','Dissolved oxygen'};

TableOutput_error=[['Output'; Output_label'],...
        ['Output Error (%)'; Output_error'],...
        ['Occurs At (h)';Occurs_at']];

save([pwd '/Results/TableOutput_error ' datestr(now)],'TableOutput_error')

end