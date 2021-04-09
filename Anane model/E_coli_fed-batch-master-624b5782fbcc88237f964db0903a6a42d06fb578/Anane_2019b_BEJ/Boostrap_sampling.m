function Y_mc = Boostrap_sampling(Pmin,data,N_samples)

%% Step-2 Check and plot residuals at the minimum fval.
% The residuals should be distributed around 0 with no apparent pattern 

[Residuals, ~, ~] = Ecoli_residual_par(Pmin,data);


[s, ~] = size(data.Offlinedata);
[t, ~] = size(data.RawDOT);

X_residuals   = Residuals(1:s);
S_residuals   = Residuals(s+1:2*s);
A_residuals   = Residuals(2*s+1:3*s);
DOT_residuals = Residuals(3*s+1:end);
% T_id = zeros(length(t_sample),1);   
%     for k1 = 1:length(t_sample) %for each sampling time
%             [~,t_idx] = min(abs(T-t_sample(k1)));
%             T_id(k1,1) = t_idx;
%     end
% DOT_residualOffl = DOT_residuals(T_id);
% Tot_res = X_residuals+S_residuals+A_residuals+DOT_residualOffl;
% All_res = vertcat(X_residuals,S_residuals,A_residuals,DOT_residualOffl);
% Plotresiduals
%% Step-2:Bootstrap sampling and re-estimation of parameters

% Solve the model using the estimated parameters. Generate synthetic data by bootstrap sampling and repeat the
% parameter estimation. New data set =  model solution + bootstrap samples

[~, Modeloutput,Time] = Ecoli_residual_par(Pmin,data);
Modeloutput_X    = Modeloutput(1:s);
Modeloutput_S    = Modeloutput(s+1:2*s);
Modeloutput_A    = Modeloutput(2*s+1:3*s);
Modeloutput_DOT  = Modeloutput(3*s+1:end);

[n,m] = size(Residuals);
NBootstrap = N_samples; 

% Pre-allocations
% Offlinedata_synX      = zeros(s,NBootstrap);
% Offlinedata_synS      = zeros(s,NBootstrap);
% Offlinedata_synA      = zeros(s,NBootstrap);
% RawDOT_syn            = zeros(length(data.RawDOT),NBootstrap);
Y_mc                  = zeros(11394,NBootstrap+1);

% Generate bootsrap samples for X, S, A and DOT separately.

for i2=1:NBootstrap  
  
    Onesample                   = ceil(s*rand(s,1));      % Generate a matrix of random (uniform) numbers, scale to datapoint range and round up. Matrix is an index used to randomly sample the residuals. 
    X_ResidualsSample           = X_residuals(Onesample);   % Random sampling with replacement from measurement errors (we respect each variable). For the relevant sample, get the residuals calculated above 
    X_SyntheticData             = Modeloutput_X+ X_ResidualsSample;  % Synthetic data or different realisations of measurement errors
%     Offlinedata_synX(:,i2)      = X_SyntheticData;
       
    S_ResidualsSample           = S_residuals(Onesample);   % Random sampling with replacement from measurement errors (we respect each variable). For the relevant sample, get the residuals calculated above 
    S_SyntheticData             = Modeloutput_S + S_ResidualsSample;  % Synthetic data or different realisations of measurement errors
%     Offlinedata_synS(:,i2)      = S_SyntheticData;
    
    A_ResidualsSample           = A_residuals(Onesample);   % Random sampling with replacement from measurement errors (we respect each variable). For the relevant sample, get the residuals calculated above 
    A_SyntheticData             = Modeloutput_A + A_ResidualsSample;  % Synthetic data or different realisations of measurement errors
%     Offlinedata_synA(:,i2)      = A_SyntheticData;

    Onesample2                  = ceil(t*rand(t,1));      % Random sampling with replacement. Generate a matrix of random (uniform) numbers, scale to datapoint range and round up. Matrix is an index used to randomly sample the residuals. 
    DOT_ResidualsSample         = DOT_residuals(Onesample2);   % Random sampling with replacement from measurement errors (we respect each variable). For the relevant sample, get the residuals calculated above 
    DOT_SyntheticData           = Modeloutput_DOT + DOT_ResidualsSample;  % Synthetic data or different realisations of measurement errors
%     RawDOT_syn(:,i2)            = DOT_SyntheticData;

    Y_mc(:,i2)             = vertcat(X_SyntheticData,S_SyntheticData,A_SyntheticData,DOT_SyntheticData);
end
    Y_mc(:,i2+1)           = Time;
    Y_mc                   =Y_mc';
end



function [Residuals,Modeloutput,Time] = Ecoli_residual_par(Pmin,data)
%Emmanuel Anane
%28/07/2016

%Simulation of fed-batch culture of E. coli. The script consists of a 
%simple batch phase, followed by an exponential fed-batch (with no product
%formation) and then a constant feed fed-batch phase.
%Load laboratory data

RawDOT = data.RawDOT;
Offlinedata = data.Offlinedata;
f       = 1;

%% Define time spans
tspan = RawDOT(:,1);
%% Define time spans
[~, FBinit]     = min(abs(tspan-11.4469));   %batch phase ends at 11.4469h, initiate exponential feed fed-batch 
[~, FBconst]    = min(abs(tspan-16.3028)); %initiate constant feed fed-batch at 16.3028
[~, Tpulse1]    = min(abs(tspan-20.5833));
[~, Tpulse2]    = min(abs(tspan-22.3639));
[~, Tpulse3]    = min(abs(tspan-24.4472));
[~, Tpulse4]    = min(abs(tspan-29.5111));
[~, Tkla]       =  min(abs(tspan-14.3417));
t_sample        = Offlinedata(:,1);
%% Initialization

y0 = [2 0.17 4.94 0.0129 98 98 0];   %Initial conditions
   %[V0, X0, S0, A0 DOT0 DOTm F0]
                                                        
Fe0 = 0.145*60/1000;             %L/h; initial value of exponential feed.

%Initial Parameter values (from literature)
      
        %[Kap     Ksa     Ko      Ks     Kia      Kis    pAmax   qAmax    qm     qSmax    Yas     Yoa     Yxa    Yem     Yos    Yxsof] 

% Par = [0.5088	0.0128	0.0001	0.0381	1.2602	1.8383	0.2286	0.1148	0.0133	0.635	0.8938	0.5221	0.5794	0.5321	1.5722	0.229]; %best Pmin
data.Par = Pmin;
%% Inputs
Si = 300;                        % concentration of glucose feed 200 g/L
mufeed = 0.222;                  % set spec. growth rate during exp. feed
DOTstar = 99;                    % equil. DO concentration at operating pressure & temperature
Kla = 220;
tau = 25;
data.u = [Si mufeed DOTstar Kla tau]; %imputs vector


%%  Batch Phase
T = [];
Y = [];
data.flag = true; %specify which equations to solve using flags
tspan1 = tspan(1:FBinit,1);
options = odeset('NonNegative',1:7);
[t1, y1] = ode15s(@fn_e_coli,tspan1,y0,options,data);
T = [T;t1];
Y = [Y;y1];
y0 = y1(end,:);

%% Exponential feed fed-batch
tspan21 = tspan(FBinit:Tkla,1);
data.flag = false; 
y0(7) = Fe0;
[t21, y21] = ode15s(@fn_e_coli,tspan21,y0,options,data);

tspan22 = tspan(Tkla:FBconst,1);
data.u(4) = 355;                             %kla increased by inc. rpm
y0 = y21(end,:);  
[t22, y22] = ode15s(@fn_e_coli,tspan22,y0,options,data);
t2 = vertcat(t21,t22(2:end));
y2 = vertcat(y21,y22(2:end,:));
%% Constant feed fed-batch, use half of the last value of exponential feed as
%   constant feed. Give intermittent pulses of glucose
data.flag = true;
y0 = y2(end,:);
y0(7) = 0.5*y0(7); 
tspan3 = tspan(FBconst:Tpulse1,1);
[t3, y3] = ode15s(@fn_e_coli,tspan3,y0,options,data);
y0 = y3(end,:);
y0(3) = 1.0; %give a pulse of 1g/L)
tspan4 = tspan(Tpulse1:Tpulse2,1);
[t4, y4] = ode15s(@fn_e_coli,tspan4,y0,options,data);
y0 = y4(end,:);
y0(3) = 0.04; %give a pulse of 0.02g/L)
tspan5 = tspan(Tpulse2:Tpulse3,1);
[t5, y5] = ode15s(@fn_e_coli,tspan5,y0,options,data);
y0 = y5(end,:);
y0(3) = 0.2; %give a pulse of 0.1g/L)
tspan6 = tspan(Tpulse3:Tpulse4,1);
[t6, y6] = ode15s(@fn_e_coli,tspan6,y0,options,data);
y0 = y6(end,:);
y0(3) = 0.2; %give a pulse of 0.08g/L)
tspan7 = tspan(Tpulse4:end,1);
[t7, y7] = ode15s(@fn_e_coli,tspan7,y0,options,data);
%% Collect all data points and extract corresponding offline data
T = [T;t2(2:end);t3(2:end);t4(2:end);t5(2:end);t6(2:end);t7(2:end)];
Y = [Y;y2(2:end,:);y3(2:end,:);y4(2:end,:);y5(2:end,:);y6(2:end,:);y7(2:end,:)];

DOT_sim     =   Y(:,6); 
c  = length(DOT_sim);
freq = 1:f:c;
DOT_sim_sampled     = DOT_sim(freq);
DOT_data            = RawDOT(:,7);
DOT_data_sampled    = DOT_data(freq);

T_id = zeros(length(t_sample),1);   
    for k3 = 1:length(t_sample) %for each sampling time
            [~, t_idx] = min(abs(T-t_sample(k3)));               
            T_id(k3,1) = t_idx;
    end
    
    X_sim = Y(T_id,2);       %locate sampling points in 
    S_sim = Y(T_id,3);       %simulated data
    A_sim = Y(T_id,4);


 Y_sim  = vertcat(X_sim, S_sim, A_sim, DOT_sim_sampled); 
 Y_data = vertcat(Offlinedata(:,2),Offlinedata(:,4),Offlinedata(:,6),DOT_data_sampled);
 
Modeloutput     = Y_sim;     
Residuals       = Y_data-Y_sim;
Time            = vertcat(t_sample,t_sample,t_sample,RawDOT(:,1));
%% Normalize the data and then the residuals
% X_simnorm           = X_sim./max(Offlinedata(:,2));
% S_simnorm           = S_sim./max(Offlinedata(:,4));
% A_simnorm           = A_sim./max(Offlinedata(:,6));
% DOT_sim_samplednorm = DOT_sim_sampled./max(DOT_data_sampled);
% 
% X_datanorm          = Offlinedata(:,2)./max(Offlinedata(:,2));
% S_datanorm          = Offlinedata(:,4)./max(Offlinedata(:,4));
% A_datanorm          = Offlinedata(:,6)./max(Offlinedata(:,6));
% DOT_data_samplednorm= DOT_data_sampled./max(DOT_data_sampled);
% 
% 
% Y_sim  = vertcat(X_simnorm, S_simnorm, A_simnorm, DOT_sim_samplednorm); 
% Y_data = vertcat(X_datanorm,S_datanorm,A_datanorm,DOT_data_samplednorm);
%  
% Modeloutput     = Y_sim;     
% Residuals       = Y_data-Y_sim;


end