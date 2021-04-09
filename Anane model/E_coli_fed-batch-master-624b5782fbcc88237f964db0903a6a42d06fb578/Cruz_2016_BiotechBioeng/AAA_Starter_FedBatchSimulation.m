
global flag1 flag2 flag3
%Emmanuel Anane
%28/07/2016

%Simulation of fed-batch culture of E. coli. The script consists of a 
%simple batch phase, followed by an exponential fed-batch (with no product
% formation) and then a constant feed fed-batch phase. Protein production is
% induced at the beginning of the constant feed fed-batch phase.
%% Define time spans
tlog = 10/3600;                      %computer logs data in 10s intervals
tspan = [0:tlog:21]';
FBinit = find(abs(tspan-14.100) < 0.001,1);  %batch phase ends at 14.1000h, initiate exponential feed fed-batch 
FBconst = find(abs(tspan-17.1667) < 0.001); %initiate constant feed fed-batch at 17.1667
%% Initialization

y0 = [2 0.0028 5 0 100 100 0 0];   %Initial conditions
   %[V0, X0, S0, A0 DOT0 DOTm F0 P0]..
                                                        
Fe0 = 0.185*60/1000;             %L/h; initial value of exponential feed.

%Initial Parameter values (from literature)
Par =[0.405 0.73 1.008 0.017 10.... % [Kap Kaq Ko Ks Ksq...
      0.46 0.65 1.3 0.04 1.357...   % [mumax pAmax qAmax qm qSmax
      1.72 0.8 0.472 1 400....      % [Yaof Yaresp Yem Yoresp Kla
      40 0.5 0.5];                  % [tau k1 k2]

%% Inputs
Si =215;            %concentration of glucose feed 200 g/L
mufeed = 0.25;       %set spec. growth rate during exp. feed
DOTstar = 99;       %equil. DO concentration at operating pressure & temperature
I = 0.45;           % mM per gram dry cell weight (inducer concentration)
Xind = 14;          %biomass at point of induction, g/L
u = [Si mufeed DOTstar I Xind]; %imputs vector
%%  Batch Phase
T = [];
Y = [];
flag1 = true; flag2 = false; flag3 = false; %specify which equations to solve using flags
tspan1 = tspan(1:FBinit,1);
options = odeset('NonNegative',[1:7]);
[t1, y1] = ode15s(@fn_e_coli,tspan1,y0,options,Par,u);
T = [T;t1];
Y = [Y;y1];
y0 = y1(end,:);
%% Exponential feed fed-batch
tspan2 = tspan(FBinit:FBconst,1);
flag1 = false; flag2 = true;
y0(7) = Fe0;
[t2, y2] = ode15s(@fn_e_coli,tspan2,y0,options,Par,u);
T = [T;t2(2:end)];
Y = [Y;y2(2:end,:)];
y0 = y2(end,:);                               
 %% Constant feed fed-batch, use the last value of exponential feed as
     %constant feed.
tspan3 = tspan(FBconst:end,1);
flag2 = false; flag3 = true;                           
[t3, y3] = ode15s(@fn_e_coli,tspan3,y0,options,Par,u);
T = [T;tspan3(2:end)];
Y = [Y;y3(2:end,:)];

%%  Plot profiles
figure
tf = {'Biomass','Substrate','Acetate','Dissolved Oxygen','Feed','Product'}; %figure titles
ly = {'X(g/L)','S(g/L)','A(g/L)','DOT(%)','F(L/h)','P(mg/L)'}; % y labels
   for i = 1:6
       w = i;
       subplot(3,2,i)
        if i >= 4
            w = i+1;
        end 
       plot(T,Y(:,w+1))
       ylabel(ly(i));  xlabel('Time (h)'); title(tf(i))
   end
   