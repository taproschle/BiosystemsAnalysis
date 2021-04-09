function [Pulse_time,V_pipette, Total_vol_added]    = Calculate_pulses(Specs,OD600,cultivation_age)   


%% Inputs for model solution
data.Par    = [0.7088	0.0128	0.0381	15.602	2.8383	0.386	0.2148	0.0133	1.235	0.8938	1.15221	0.5794	0.5321	1.2722	0.229  0.35]; %best parameter values from model already fitted to strain

DOTstar     = 100;
Kla         = 800;
tau         = 15;
gamma       = 0.8;                             
F           = 0;
data.u      = [Specs.Si  DOTstar Kla tau gamma F 0];
data.flag   = 3;
%% Define pulse feed parameters
number_reactors = length(OD600);
mu_set          = Specs.mu_set;
X               = abs(OD600.*0.37);                     % Biomass measured as OD600, convert to CDW. E coli BW25113 has conversion factor of 0.37
V               = Specs.V_reactor./1000;                % convert reactor volumes to L
tlog            = 10/3600;
feed_time       = Specs.T_vector(1);                 
resting_time    = Specs.T_vector(2);                
T_start         = cultivation_age;       

cycle_time      = feed_time+resting_time;

deltse          = cycle_time/60;               
                                      
deltsp          = feed_time/60;                
                               
Fe0             = X.*V.*(mu_set./(data.Par(13)*Specs.Si));        % Calculate initial feed rate to maintain mu_set, L/h

%% Exponential Pulse fed-batch
Fpulse          = zeros(Specs.pulses,number_reactors);            % flow rate in case of pulse feed, L/h
V_pulse         = zeros(Specs.pulses,number_reactors);            % volume of pulses, what the pipette draws per feed cycle
X_next          = zeros(1,number_reactors);
Volumes         = zeros(Specs.pulses,number_reactors);            % Bioreactor volume in ml
Pulse_time      = zeros(Specs.pulses,1); 
w0              = zeros(number_reactors,6);                 % 6 is the number of state variables in y0
w0(:,[1;4;5])   = [X',100*ones(number_reactors,2)];         % Initial conditions
V_reactor       = Specs.V_reactor;


for k1  = 1:Specs.exponential_time*60/sum(Specs.T_vector);                                             % first 3hrs of pulses are done under exponential feeding.  
    r = k1-Specs.Sample_times;                                    % determine sampling times
    if any(r==0)
        Sample = 1;                                    
    else
        Sample = 0;
    end
    
    tse = [T_start T_start+deltse];                     %integrate Fe in this interval(10min).
    tsp = tse(1):tlog:tse(1)+deltsp;                    %Feed equival.Vpulse in this interval(1min).
    fe = @(t) Fe0*exp(mu_set*t);                        %Function of the exponential feed rate, fe. 
    Intfe = integral(fe,0,deltse,'ArrayValued',true);   %integrate fe to find area under the exponential feed profile. This is equal to the total volume of feed given in 10min. The integral used here makes the calculations slower, compared to calculating the analytical solution of fe. It is maintained just for clarity, since the code is solved only once, per fermentation.
    V_pulse(k1,:) = max(5,round(Intfe.*(1000*1000)));                 %give this volume in 1 minute, convert to uL.             % L
    Fpulse(k1,:)  = Intfe/deltsp;                       % calculate the corresponding feed flow rate in L/h, used only for the model solution.

  % Solve the model to get the new biomass to be used to calculate the pulse size, reconcile volumes and calculate dilution rates and new concentrations.
    
    for k2 = 1:number_reactors
   
        y0            = w0(k2,:);
        data.u(6)     = Fpulse(k1,k2);
        data.u(7)     = V(1,k2);                        % select which reactor to use
        options       = odeset('NonNegative',1:6);
        [~, y2]       = ode15s(@fn_Ecoli,tsp,y0,options,data);

        y02           = y2(end,:);

        data.u(6)     = 0;                              % feed is switched off for 9 min
        tswof         = tsp(end):tlog:tse(end);         % tswof = time span without feed(9 min).
        [~, y3]       = ode15s(@fn_Ecoli,tswof,y02,options,data);

        w0(k2,:)      = y3(end,:);                      % set new parameters for looping.

        X_next(1,k2)  = y3(end,1);                      % new biomass concentration.

        V_reactor(1,k2) = V_reactor(1,k2) + V_pulse(k1,k2)./1e3 - Sample.*Specs.V_sample(1,k2);     %% Calculate the new volume after pulse addition and sampling. Sample is either 1 (for sampling) or 0 (no sampling), V_sample contains the sampling volumes of all 8 bioreacors

    end
    V                  = V_reactor./1000;
    Volumes(k1,:)      = V;
    X_new              = w0(:,1)';
    Fe0                = X_new.*V.*(mu_set./(data.Par(13)*Specs.Si));        % L/h
    T_start            = tse(end);
    Pulse_time(k1)    = tse(1);
end


%% Constant pulse, using the last value of the exponential feed as the constant pulse value

    for k3  = k1+1:Specs.pulses

        r   = k3-Specs.Sample_times;                               % determine time samples are taken
            if any(r==0)
                Sample = 1;
              else
                Sample = 0;
            end

        Fpulse(k3,:)  = Intfe./deltsp;                                   
        V_pulse(k3,:) = max(5,round(Intfe.*(1000*1000)));    %minimum volume that the Tecan can pipette is 5 micro liters, if the calculated value is less, the 5uL is taken as V_pipette                 
        for k4 = 1:number_reactors
            V_reactor(1,k4) = V_reactor(1,k4) + V_pulse(k3,k4)./1e3 - Sample.*Specs.V_sample(1,k4);        % Sample is either 1 (for sampling) or 0 (no sampling), V_sample contains the sampling volumes of all 8 bioreacors
            
        end
        V                   = V_reactor./1000;
        Volumes(k3,:)       = V; 
        Pulse_time(k3)     = Pulse_time(k3-1)+deltse;
        
%             if Sample == 1;
%                 volumes_ratio = Volumes(k3,:)./Volumes(k3-1,:);
%                 Intfe = Intfe.*volumes_ratio;
%             end

    end
    

    Pulse_time             = Pulse_time.*3600;

    V_pipette               = V_pulse;                           % Volume of feed to be pipetted (in microliters)


%% Calculate pulse size for bioreactors with double (or multiple) intervals between pulses
    
    if strcmp(Specs.double_pulses_int,'yes');

     Volume_db_pulse = double_pulses(Specs,OD600,cultivation_age);  
     
     [p, q] = size(Volume_db_pulse);
     V_pipette(:,Specs.which_bioreactors_dbp) = 0;
     
         for k4 = 1:q 
                  rxtor = Specs.which_bioreactors_dbp(k4);
             for k5 = 1:p 
                 V_pipette(k5*Specs.Pulse_interval,rxtor) = round(Volume_db_pulse(k5,k4));
             end
         end
    end
    
%% Remove setpoints for bioreactors with no pulses
    
    if strcmp(Specs.no_pulses,'yes');

     V_pipette(:,Specs.which_bioreactors_np) = 0;
     
    end

%% Add enzyme solution in pulses
    
    if strcmp(Specs.enzyme_addition,'yes');

     
     V_pipette(:,Specs.which_bioreactors_enz) = max(5,round(V_pipette(:,Specs.which_bioreactors_enz)./Specs.enz_Vfactor));
     
         
    end


%% Plot feed profiles
    
    Total_vol_added         = sum(V_pipette)./1000;

    if strcmp(Specs.plot_feed_profile,'yes');
    figure
%     subplot(1,2,1)
    plot(V_pipette(:,:),'o')
    ylabel('Pulse volume (uL)'); xlabel('Number of Pulses x 5 (min)')
    set(gca,'LineWidth',2,'FontSize',14,'FontWeight','bold') 
   
    figure
%     subplot(1,2,1)
    bar(V_pipette)
    ylabel('Pulse volume (uL)'); xlabel('Number of Pulses x 5 (min)')
    set(gca,'LineWidth',2,'FontSize',14,'FontWeight','bold') 

%     subplot(1,2,2)
    figure
    plot(1000.*Volumes(:,:))
    ylabel('Bioreactor Volume (mL)'); xlabel('Number of Pulses x 5 (min)')
    set(gca,'LineWidth',2,'FontSize',14,'FontWeight','bold') 
    end   
end

function dydt= fn_Ecoli(t, y, data)
%{ 
Emmanuel Anane 
Institute fur Biotechnologie, TU-Berlin.
09/11/2016

this function is the set of ODEs for a fed-batch culture of e coli. y(1) = dV/dt,
y(2) = dx/dt, y(3) = dS/dt, y(4) = dP/dt, y(5) = dF/dt. The script that 
calls it implements a fed-batch process of e. coli fermentation, with
product formation after iptg induction.
%}
%% 

if isempty(t)       %for calculating Jacobians
    
    syms('V', 'X', 'S', 'A', 'DOT', 'DOTm', 'F', 'P','Cs', 'Cx', 'H','gamma',...
        'Si','DOTstar','Kla','tau')
    
    Cs          = 0.391;                   
    Cx          = 0.488;                     
    H           = 14000;    
    Ko          = 0.001;
else
    %% State variables
    X       = y(1);           %biomass
    S       = y(2);           %substrate
    A       = y(3);           %acetate
    DOTa    = y(4);           %dissolved oxygen, no probe response
    DOTm    = y(5);           %dissolved oxygen, with probe response time         
    P       = y(6);
    %%  Fermentation constants
    Cs  = 0.391;               %carbon content of glucose (ref:enfors)
    Cx  = 0.488;               %carbon content of biomas(e coli)  
    H   = 14000;               %Henry's law constant to convert from mol to % DO
    Ko  = 0.001;

    %% INPUTS
    Si          = data.u(1);        %feed concentration (g/L)
    DOTstar     = data.u(2);     
    Kla         = data.u(3);
    tau         = data.u(4);
%     gamma       = data.u(5);
    F           = data.u(6);
    V           = data.u(7);
end
%%  Parameters

if isempty(t)
    syms('Kap', 'Ksa', 'Ko', 'Ks', 'Kia', 'Kis', 'pAmax', 'qAmax', 'qm',...
        'qm', 'qSmax', 'Yas', 'Yoa', 'Yxa', 'Yem', 'Yos', 'Yxsof', 'Ypx','k')
else
    
    Kap     = data.Par(1);   %monod-type saturation constant, intracellular acetate production
    Ksa     = data.Par(2);   %affinity constant, acetate consumption (0.05 g/L)
    Ks      = data.Par(3);   %affinity constant, substrate consumption (0.05 gglu./L)
    Kia     = data.Par(4);   %inhibition constant, inhib. of glucose uptake by acetate
    Kis     = data.Par(5);   %inhibition constant, inhib. of ace.uptake by glucose(4g/L)
    pAmax   = data.Par(6);   %max spec acetate production rate (0.15 g_ace/(gx.h)))
    qAmax   = data.Par(7);   %max spec acetate consumption rate (0.15 g_ace/(gx.h)))
    qm      = data.Par(8);   %spec maintenance coefficient (0.04 g_glu/(gx.h))
    qSmax   = data.Par(9);   %max spec glucose uptake rate (1.3 g_glu./gx.h)
    Yas     = data.Par(10);  %yield of acetate on substrate (g.ace/g.glu)
    Yoa     = data.Par(11);  %yield of oxygen on acetate g/g)
    Yxa     = data.Par(12);  %yield of biomass on acetate
    Yem     = data.Par(13);  %yield exclusive maintenance (Yem = Yxs-qm)
    Yos     = data.Par(14);  %yield of oxygen on glucose (g/g)
    Yxsof   = data.Par(15);  %yield of biomass on all other products of overflow routes, excluding acetate
    Ypx     = data.Par (16); %yield of product on biomass (g./g.h)
%     k       = data.Par(18);  % non-growth associated product formation rate (g/g.h)
end
%% Explicit algebraic equations
qS      = (qSmax./(1+(A/Kia)))*(S/(S+Ks));
qSof    = pAmax*(qS/(qS+Kap));
pA      = qSof*Yas;                             %acetate production from the overflow flux, and not from total qS
qSox    = (qS-qSof)*(DOTa/(DOTa+Ko));             %When there is no oxygen, only the qsof is operational            
qSan    = (qSox-qm)*Yem*Cx/Cs;      
qsA     = (qAmax./(1+(qS/Kis)))*(A/(A+Ksa));    %acetate consumption, what happens if A is replaced with pA?
qA      = pA-qsA;                               %net acetate accumulation
mu      = (qSox-qm)*Yem + qsA*Yxa + (qSof-pA)*Yxsof;
qO      = Yos*(qSox-qSan)+qsA*Yoa;
Kp      = (1/tau)*3600;
% qp      = ((mu*Ypx) + k).*gamma;
qp      = mu.*Ypx;%.*exp(gamma);

DOT     = (Kla*DOTstar-qO*X*H)./Kla;

%% ODE system
dXdt    = X*(mu - (F/V));             %dX/dt
dSdt    = F/V*(Si - S) - (qS*X);      %dS/dt
dAdt    = qA*X-(F/V)*A;               %dA/dt
dDOTadt = Kla*(DOTstar-DOTa)-qO*X*H;   %dDOT/dt
dDOTmdt = Kp*(DOT-DOTm);              %dissolved oxygen with probe response 

    if isempty(t);                       %use this to generate Jacs
           if data.flag == 1              %batch phase
              dPdt = 0;
              elseif data.flag == 2     %exponential fed-batch phase, non-production.
              dPdt = 0;
              elseif data.flag == 3;    %constant feed fed-batch phase
              dPdt = qp - mu*P;
              
            end
    else  
           if data.flag == 1                       
                dPdt = 0;
                elseif data.flag == 2   
                dPdt = 0;
                elseif data.flag == 3;  
                dPdt = qp - mu*P;
            end
    end
   if ~isempty(t)
        if length(y) == 6 %
            dydt = [dXdt;dSdt;dAdt;dDOTadt;dDOTmdt;dPdt];
        else
            m = length(JacX); % number of state variables
            n = length(JacP);%./length(JacX); % number of Parameters
            Sr = reshape(y(7:end),[m,n]); 
            Sm = JacX*Sr;       %relative sensitivity matrix
            SLi = reshape(Sm,m*n,1);
            JacPe = reshape(JacP,[m*n,1]);
            dydt = [dXdt;dSdt;dAdt;dDOTadt;dDOTmdt;dPdt;JacPe+SLi];
        end
    else
      dydt = [dXdt;dSdt;dAdt;dDOTadt;dDOTmdt;dPdt];
    end
end
