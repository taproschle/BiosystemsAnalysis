function Z = Objective_function(Par_norm,data)



Offlinedata_vec     = data.Offlinedata_vec;
RawDOT              = data.RawDOT;
freq                = data.freq;
wt                  = data.wt;
off_idx             = data.off_idx;
Fe0                 = data.Fe0;
%% Define time spans
tspan        = RawDOT(:,1);
[~, FBinit]     = min(abs(tspan-11.4469));        % batch phase ends at 11.4469h, initiate exponential feed fed-batch 
[~, FBconst]    = min(abs(tspan-16.3028));       % initiate constant feed fed-batch at 16.3028
[~, Tpulse1]    = min(abs(tspan-20.5833));
[~, Tpulse2]    = min(abs(tspan-22.3639));
[~, Tpulse3]    = min(abs(tspan-24.4472));
[~, Tpulse4]    = min(abs(tspan-29.5111));
[~, Tkla]       = min(abs(tspan-14.3417));
t_sample        = Offlinedata_vec(:,1);
%% Initialization
y0          = data.y0;
data.Par    = Par_norm.*data.Par_ref;
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

%% Visualize
% Plotprofiles
%% Extract and prepare data from offline measurments
[a, ~] = size(Y);               %for the length of the simulated data
z_points = 1:freq:a;            %set frequency at which obj function points should be taken.


X_data = Offlinedata_vec(off_idx(1):off_idx(2),2);
X_normscaled = wt(1).*(X_data./(max(X_data)));   %normalize data to [0 1]

S_data = Offlinedata_vec(off_idx(3):off_idx(4),2);  
S_normscaled = wt(2).*(S_data./(max(S_data)));

A_data = Offlinedata_vec(off_idx(5):off_idx(6),2);
A_normscaled = wt(3).*(A_data./(max(A_data)));

DOT_data = RawDOT(:,2);
DOT_data = DOT_data(z_points);
DOT_normscaled = wt(4).*(DOT_data./(max(DOT_data)));


Mess_data = vertcat(X_normscaled,S_normscaled,A_normscaled,...
                    DOT_normscaled);

%% Extract corresponding simulations

T_id = zeros(length(t_sample),1);   
    for k3 = 1:length(t_sample) %for each sampling time
            [~, t_idx] = min(abs(T-t_sample(k3)));               
            T_id(k3,1) = t_idx;
    end
 
X_sim = Y(T_id(off_idx(1):off_idx(2)),2);           %locate sampling points in
X_sim_normSc = wt(1).*(X_sim./(max(X_data)));          %normalized and scaled data

S_sim = Y(T_id(off_idx(3):off_idx(4)),3);           %locate sampling points in
S_sim_normSc = wt(1).*(S_sim./(max(S_data)));          %normalized and scaled data

A_sim = Y(T_id(off_idx(3):off_idx(4)),4);
A_sim_normSc = wt(3).*(A_sim./(max(A_data)));

DOT_sim  = Y(:,6);
DOT_sim  = DOT_sim(z_points);   %select points for calculating objective function.
DOT_sim_normSc  = wt(4).*(DOT_sim./(max(DOT_data))); 


Sim = vertcat(X_sim_normSc,S_sim_normSc,A_sim_normSc,...
              DOT_sim_normSc);
%% Regularization (Tikhonov)
if(data.Tikh_Reg)
    Tikh = data.Par-data.Par_ref;
    Reg_val = data.Tikh_lambda.*sum(Tikh.^2);
else
    Reg_val = 0;

end

    %% Define objective function (least square problem) for PE
  if data.optroute == 0
         try
                    Z = sum((Mess_data - Sim)'*(Mess_data - Sim))+ Reg_val;            %for fmincon
         catch                
                    Z = 1e5.*sum(Mess_data);
                    disp('Solver failed, using dummy Z')
         end
   else
          try
                Z = Mess_data - Sim;            %for lsqnonlin
        catch                
                Z = 1e5.*Mess_data;
                disp('Solver failed, using dummy Z')
         end
  end
end