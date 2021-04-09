%{ 

%Emmanuel Anane

Using CEDEX and LABCHIP  data d
 
23 May 2018

Description: Simulation of fed-batch culture of E. coli. The script 
implements a pulse feeding profile after an initial batch phase as follows: 
(i) FBexp--Exponential fed-batch where the feed is divided into discrete
pulses, and a pulse is injected into the bioreactors every 5 minutes. 
During feeding, the mass of glucose delivered is a pulse is equal to what 
would have been fed exponentially in a continuous feed for the next 5 minutes 
when the feed is off. This creates glucose concentration gradients. 
(ii) FBconst--the last value of the exponential pulse is maintained till 
the end of tspan. Pulse timing is the same as for FBexp, except that the 
pulse size remains constant. Protein production is started by IPRG
addition to a final concentration of 0.5mM.


Kla change at 8.3h within exponential fed-batch phase.
%}
%% move to PE

% dbstop if error
load mbp6.mat; load Pmin.mat ;load Pulse_time_2XP.mat
Rxtor_num           = 14;
tspan               = Bioreactor(Rxtor_num).DOT(25:5040,1);
DOT                 = Bioreactor(Rxtor_num).DOT(25:5040,2);

% Rx_num = [3 11 19]; %1xP cultivation
Rx_num = [6 14 22]; %2xP cultivation
% Rx_num = [13 21]; %1xP Enbase cultivation

data.Offlinedata = [];

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
     
[~, FBinit]                = min(abs(tspan-6.2956));  %these time points are valid for all 24 bioreactors 

kla_change(1)              = 1;   

[~, kla_change(2)]         = min(abs(tspan-2.957));   

[~, kla_change(3)]         = min(abs(tspan-3.819));   

[~, kla_change(4)]         = min(abs(tspan-4.448));   

[~, kla_change(5)]         = min(abs(tspan-4.943));   

kla_change(6)              = FBinit;   

% t_sample            = Offlinedata1(:,1);
%% Inputs

Si                  = 200;
F                   = 0;
V                   = 0.01;
mu_set              = 0.35;

data.u              = [Si  V F mu_set];
   
%% Initialization

for b1 = 1:1
% y0               = data.y0;
% y0          = [0.037 5   0  DOT(1)  0];


                %[Kap     pAmax   qAmax    qSmax   Yxa     Yem     Yxsof   Yoa     Yos] 
% data.Par     = [0.860  0.01354 0.00487 0.878 1.5802   0.3107  0.50737 0.04094  1.9489 0.24566  0.22186  0.44906 0.19788 0.54244  0.95752 1.16299  350  160  220  300 420  100]; 
% data.Par    =  [1.0587  0.03805  0.0598 1.5197   3.62811   0.7843   0.1701   0.0607  1.937  0.5761    0.5304   0.5121  0.4134   0.3934  1.2718  1.276 643   573   701  432  346]; % PE1
% data.Par      = [0.1169  0.7074  0.2834  1.634 0.6893 0.4155 0.56136  1.217 1.441]; %PE2

              %Ksa       Ks         Kip    Kis    qm     Yas     Ypx     Kla
% data.Konst  = [0.0603    0.07588   5.646  8.869  0.06635  0.7803  0.5773 476];
                
                %[Kap     Ksa     Ks    Kip    Kis     pAmax   qAmax  qm   qSmax   Yas     Yxa     Yem     Yxsof   Ypx  Yoa     Yos tau kla] 

% data.Par      = [0.1169  0.0603 0.07588 5.646  8.869  0.7074  0.2834  0.06635 1.634 0.9803  0.1893 0.5155 0.16136  0.5773 1.217 1.441]; %PE2
% data.Par       = [0.2580 0.0498 0.0563 4.3444 6.3751 0.5898 0.1033 0.0421 1.3066 0.5356 0.1503 0.5362 0.8843 0.3157 2.3628 2.4953 10  475  523]; 
% data.Par = [0.497412727793123,0.0200617814679632,0.0398372570756547,3.95678736196967,1.30414757914609,0.328112624523344,0.117964671922563,0.0337145574372697,2.32914847523908,0.732498809511869,0.212408160528296,0.671372358963235,0.216756510439669,0.349664078949334,2.76349978014211,2.35002387046722,47.9376562590565,433.979797008276,530.537479823478];
% data.Par = [0.497412727793123,0.0200617814679632,0.0398372570756547,3.95678736196967,1.30414757914609,0.428112624523344,0.117964671922563,0.0337145574372697,2.32914847523908,0.732498809511869,0.212408160528296,0.671372358963235,0.216756510439669,0.349664078949334,2.76349978014211,2.35002387046722,47.9376562590565,433.979797008276,530.537479823478];

data.Par = Pmin(b1,:);
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
y0          = [2.804 0   0.105  DOT(FBinit)  0];

N_pulses    = length(Pulse_time);
Exp_pulses  = 15;
data.flag   = 1; 
Fe0         = 9.065e-5;  %from experiment
deltsp      = 1/60;
deltse      = 10/60;             %Time span to calculate the exp. feed.   
Pulse = Pulse_time';
data.change = 1;
for k1          = 1:Exp_pulses
    
        if k1>12
            data.change = 2; % Kla change after the 24th pulse
        end
    
    T_start     = Pulse(1,k1);               %Time to start the pulse
    mu_set      = data.u(4);
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
y0(5)       = Bioreactor(Rxtor_num).Product(2,2)./1000;

for k2 = 1:N_pulses-Exp_pulses-1;
    
    w       = k1+k2;  

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

Plotprofiles

end



