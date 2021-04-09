function Z = ObjFn_mbp6(Par_norm,data)


%% Initialization

% freq            = data.freq;
wt              = data.wt;
off_idx         = data.off_idx;

y0              = data.y0;
              
data.Par        = Par_norm.*data.Par_ref;

tspan           = data.tspan;
kla_change      = data.kla_change;
Offlinedata     = data.Offlinedata;
t_sample        = Offlinedata(:,1);
DOT             = data.DOT;
Pulse_time      = data.Pulse_time;
%%  Batch Phase
T = [];
Y = [];
data.flag   = 1;               %specify which equations to solve using flags
options     = odeset('NonNegative',1:5);

% for k0 = 1:5; %Kla changed 5 times in cultivation
% data.change = k0;
% tspan1      = tspan(kla_change(k0):kla_change(k0+1));
% [t1, y1]    = ode15s(@fn_Ecoli,tspan1,y0,options,data);
% 
%     if k0 == 1
%         
%         T           = [T;t1];
%         Y           = [Y;y1];
%     else
%         T           = [T;t1(2:end)];
%         Y           = [Y;y1(2:end,:)];
%     end
%     
%     y0          = y1(end,:);
% 
% end
% plot(T,Y(:,4));hold on
% plot(tspan(2:FBinit),DOT(2:FBinit));hold off


%% Exponential Pulse fed-batch 

N_pulses    = length(Pulse_time);
Exp_pulses  = 30;
data.flag   = 1; 
Fe0          = 9.065e-5;  %from experiment
deltsp       = 1/60;
deltse      = 5/60;             %Time span to calculate the exp. feed.   
mu_set      = data.u(4);
V           = data.u(2);
Pulse = Pulse_time';
data.change = 1;
for k1          = 1:Exp_pulses
    
        if k1>13 && k1<=24
            data.change = 2; % Kla change after the 24th pulse
        elseif k1>24 && k1<=34
            data.change = 3; % Kla change after the 34th pulse 
       end
        
    
    T_start     = Pulse(1,k1);               %Time to start the pulse
    [~, p_start]         = min(abs(tspan-T_start));      
    [~, p_end]           = min(abs(tspan-Pulse(1,k1+1)));   

    tsp = tspan(p_start:p_end); %solve ode system for this pulse time.
    fe      = Fe0*exp(mu_set*deltse);         %Analytical solution of F. 
    V_pulse   = (1./mu_set).*Fe0.*exp(mu_set*deltse)-...
    (1./mu_set)*Fe0*exp(mu_set*0);       %integrate to find area under F.
    fp      = V_pulse./deltsp;                  %define fp as area/time (L/h)
                                         %F = Par.Fp in the function to 
                                         %include Fp in the ODE parameters.

    data.u(3)    = fp;                                       % set pulse feed rate;
    Glu_in_pulse = (V_pulse*data.u(1) + y0(2)*V)./(V_pulse+V);
    y0(2) = Glu_in_pulse;
    
    V            = V + V_pulse;                    % re-calculate volume.
    data.u(2)    = V;

    options      = odeset('NonNegative',1:5);
   [t2, y2]      = ode15s(@fn_Ecoli,tsp,y0,options,data);

    y0           = y2(end,:);
    if k1 == 1     
        T           = [T;t2];      %concatenate t's and y's.
        Y           = [Y;y2];
    else
        T           = [T;t2(2:end)];      %concatenate t's and y's.
        Y           = [Y;y2(2:end,:)];
    end
    Fe0 = fe;
end
%% Constant pulse fed-batch
data.flag   = 2; 
y0(5)  = Offlinedata(off_idx(7),2)/1000; % product concentration at the time of induction
for k2 = 1:N_pulses-Exp_pulses-1;
    
    w       = k1+k2;  

        if w>34 
            data.change = 4; % Kla change after the 24th pulse
        end

    T_start = Pulse(1,w);
    
    [~, p_start]         = min(abs(tspan-T_start));      
    [~, p_end]           = min(abs(tspan-Pulse(1,w+1)));   

    tsp = tspan(p_start:p_end); %solve ode system for this pulse time.

    
    data.u(3)   = fp;                               % set constant pulse feed rate;
    y0(2) = Glu_in_pulse;


    V            = V + V_pulse;                    % re-calculate volume.
    data.u(2)   = V;

    [t4, y4]    = ode15s(@fn_Ecoli,tsp,y0,options,data);

    y0          = y4(end,:);
    
    T           = [T;t4(2:end)];         
    Y           = [Y;y4(2:end,:)];
end

%% Visualize profiles
   if strcmp(data.PE_status,'no')

     Plotprofiles
     Z = 0;

     if strcmp(data.Metabolic_states,'yes')
         fn_celstate(T,Y,data)
     end
   else

    %% Extract and prepare data from offline measurments
    % [a, ~]      = size(Y);               %for the length of the simulated data
    % z_points    = 1:freq:a;            %set frequency at which obj function points should be taken.

    X_data      = 0.37.*Offlinedata(off_idx(1):off_idx(2),2); %convert ODs to dry weight
    X_normscaled = wt(1).*(X_data./(max(X_data)));   %normalize data to [0 1]

    S_data = Offlinedata(off_idx(3):off_idx(4),2);  
    S_normscaled = wt(2).*(S_data./(max(S_data)));

    A_data = Offlinedata(off_idx(5):off_idx(6),2);
    A_normscaled = wt(3).*(A_data./(max(A_data)));

    DOT_data = DOT(kla_change(6):end);
    % DOT_data = DOT_data(z_points);
    DOT_normscaled = wt(4).*(DOT_data./(max(DOT_data)));

    P_data = Offlinedata(off_idx(7):off_idx(8),2)./1000;
    P_normscaled = wt(5).*(P_data./(max(P_data)));   %normalize data to [0 1]



    Mess_data = vertcat(X_normscaled,S_normscaled,A_normscaled,...
                        DOT_normscaled,P_normscaled);

    %% Extract corresponding simulations

    T_id = zeros(length(t_sample),1);   
        for k3 = 1:length(t_sample) %for each sampling time
                [~, t_idx] = min(abs(T-t_sample(k3)));               
                T_id(k3,1) = t_idx;
        end

    X_sim = Y(T_id(off_idx(1):off_idx(2)),1);           %locate sampling points
    X_sim_normSc = wt(1).*(X_sim./(max(X_data)));          %normalized and scaled data

    S_sim = Y(T_id(off_idx(3):off_idx(4)),2);           %locate sampling points
    S_sim_normSc = wt(2).*(S_sim./(max(S_data)));          %normalized and scaled data

    A_sim = Y(T_id(off_idx(5):off_idx(6)),3);
    A_sim_normSc = wt(3).*(A_sim./(max(A_data)));

    DOT_sim  = Y(:,4);
    % DOT_sim  = DOT_sim(z_points);   %select points for calculating objective function.
    DOT_sim_normSc  = wt(4).*(DOT_sim./(max(DOT_data))); 

    P_sim = Y(T_id(off_idx(7):off_idx(8)),5);      
    P_sim_normSc = wt(5).*(P_sim./(max(P_data)));


    Sim = vertcat(X_sim_normSc,S_sim_normSc,A_sim_normSc,...
                  DOT_sim_normSc, P_sim_normSc);

        %% Define objective function (least square problem) for PE
      if data.optroute == 0
             try
                        Z = sum((Mess_data - Sim)'*(Mess_data - Sim));    %for fmincon
             catch                
                        Z = 1e2.*sum(Mess_data);
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
end