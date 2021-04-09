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
the end of tspanOrig. Pulse timing is the same as for FBexp, except that the 
pulse size remains constant. Protein production is started by IPRG
addition to a final concentration of 0.5mM.
%}
%% move to PE

% dbstop if error
load mbp6.mat ; load Pmin.mat; load Pulse_time_12

Rxtor_num           = 10;
tspanOrig               = Bioreactor(Rxtor_num).DOT(25:4970,1);
DOT                 = Bioreactor(Rxtor_num).DOT(25:4970,2);

Rx_num = [2 10 18]; %triplicate reference cultivation

 %% Extract event times during cultivation                            
    
[~, FBinit]                = min(abs(tspanOrig-5.30122));  %these time points are valid for all 24 bioreactors 

kla_change(1)              = 1;   

[~, kla_change(2)]         = min(abs(tspanOrig-2.957));   

[~, kla_change(3)]         = min(abs(tspanOrig-3.819));   

[~, kla_change(4)]         = min(abs(tspanOrig-4.448));   

[~, kla_change(5)]         = min(abs(tspanOrig-4.943));   

kla_change(6)              = FBinit;

[~, kla_change(7)]         = min(abs(tspanOrig-8.3036));

[~, TInduce]               = min(abs(tspanOrig-8.93333));  %these time points are valid for all 24 bioreactors 

% t_sample            = Offlinedata1(:,1);

%% for fed-batch phase only
tspan               = tspanOrig(FBinit:end);
DOT                 = DOT(FBinit:end);




%% collect data of all triplicate bioreactors of this condition
data.Offlinedata = [];

for k0 = 1:length(Rx_num)
    
Plot_data(k0).biomass =  Bioreactor(Rx_num(k0)).OD600(2:end,:);

data.Offlinedata         = [data.Offlinedata;...
                           Plot_data(k0).biomass];
                    
end

for k0 = 1:length(Rx_num)
Plot_data(k0).Glu =  Bioreactor(Rx_num(k0)).Glucose(2:end,:);

data.Offlinedata        = [data.Offlinedata;...
                            Plot_data(k0).Glu];
end
                            
for k0 = 1:length(Rx_num)
Plot_data(k0).Ace =  Bioreactor(Rx_num(k0)).Acetate(2:end,:);

data.Offlinedata        = [data.Offlinedata;...
                            Plot_data(k0).Ace];
end

for k0 = 1:length(Rx_num)
Plot_data(k0).Pro =  Bioreactor(Rx_num(k0)).Product(2:end,:);
    
data.Offlinedata        = [data.Offlinedata;...
                            Plot_data(k0).Pro];
end                            
    
%% Inputs

Si                  = 200;
F                   = 0;
V                   = 0.01;
mu_set              = 0.35;
pH                  = 6.87;
U0                  = 0.2;
data.u              = [Si  V F mu_set pH];
   
%% Initialization
for b1 = 1:1

% y0               = data.y0;
y0          = [0.037 5   0  DOT(1)  0  0];

            % [Ksa     Ks    Kip   Kis    qm     Yas    kla1 kla2 kla3 kla4 kla5
data.Konst = [0.0472 0.0334 2.9640 7.875  0.0202 0.9017];

% data.Par   = [0.735953118517400,0.615927843922152,0.329795153683726,2.08430470379284,0.661943697496441,0.422566917299501,0.430882092629384,0.424285612611608,1.01737938814766,1.32707658179553,1.12726870623932,6.82196122263109 474 342  437  463 651];

            %[Kap     pAmax   qAmax    qSmax     Yxa     Yem     Yxsof    Ypx     Yoa     Yos      U1      U2] 
            
% data.Par     = [0.860  0.01354 0.00487 0.878 1.5802   0.3107  0.50737 0.04094  1.9489 0.24566  0.22186  0.44906 0.19788 0.54244  0.95752 1.16299  1.30   7.5  300  200  300  350 400  500];

% data.Par = [0.6428 0.6453 0.3299 2.0468 0.7524 0.4342 0.4304 0.3918 0.9977 1.3801 1.0385 6.4228 515 294 369 482 589 780];

% data.Par = [1.05561059766142,0.584191027013514,0.314835585139339,2.00737020473958,0.735743266229447,0.446672273177751,0.406217795212519,0.424148449093979,1.10697681254694,1.43169319807844,1.02202593939202,6.14675831747499 564 780];
data.Par = Pmin(b1,:);
%%  Batch Phase
T = [];
Y = [];
% data.flag   = 0;               %specify which equations to solve using flags
% options     = odeset('NonNegative',1:5);
% 
% for k0 = 1:5; %Kla changed 5 times in cultivation
% data.change = k0;
% tspanOrig1      = tspanOrig(kla_change(k0):kla_change(k0+1));
% [t1, y1]    = ode15s(@fn_Ecoli,tspanOrig1,y0,options,data);
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
% 
% plot(T,Y(:,4));hold on
% plot(tspanOrig(2:FBinit),DOT(2:FBinit));hold off

%% Fed-batch phase: non-induced 
data.flag = 1;
data.change = 5;
y0 = [2.56,0,0.49,93.17,0,0];
y0(6)        = U0; % initial glucose release rate
tspanOrig21      = tspanOrig(FBinit:kla_change(7));
options      = odeset('NonNegative',1:5);
[t21, y21]   = ode15s(@fn_Ecoli,tspanOrig21,y0,options,data);
   
y0          = y21(end,:);

data.change = 6;
y0(6) = 0.6;
tspanOrig22      = tspanOrig(kla_change(7):TInduce);

[t22, y22]    = ode15s(@fn_Ecoli,tspanOrig22,y0,options,data);
   
y0          = y22(end,:);


T           = [T;t21(1:end);t22(2:end)];      %concatenate t's and y's.
Y           = [Y;y21(1:end,:);y22(2:end,:)]; 

% plot(T,Y(:,4));hold on
% plot(tspanOrig(2:TInduce),DOT(2:TInduce));hold off

%% Constant pulse fed-batch: protein production phase
data.flag   = 2; 

tspanOrig3      = tspanOrig(TInduce:end);

[t3, y3]    = ode15s(@fn_Ecoli,tspanOrig3,y0,options,data);


T           = [T;t3(2:end)];         
Y           = [Y;y3(2:end,:)];
    %% Visualize profiles

Plotprofiles
end




