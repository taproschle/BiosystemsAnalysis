
global flag1 flag2 flag3

%Emmanuel Anane, Diana C Lopez C, Peter Neubauer, M Nicolas Cruz Bournazou
%Chair of Bioprocess Enginering,
%Institute fur Biotechnologie, TU-Berlin.
%28/01/2017


% Model for the paper 'Modelling overflow metabolism in E. coli by acetate
% cycling', published in Biochemical Engineering Journal...(Vol, page).


%Simulation of fed-batch culture of E. coli. The script consists of a 
%simple batch phase, followed by an exponential fed-batch (with no product
%formation) and then a constant feed fed-batch phase. Intermittent glucose
%pulses are given during the constant feed to capture the dynamics of
%glucose consumption and its conversion to acetate in overflow metabolism.

%% Load laboratory data

load RawDOT
load Offlinedata


%% Define time spans

tspan        = RawDOT(:,1);
[~, FBinit]  = min(abs(tspan-11.4469));   %batch phase ends at 11.4469h, initiate exponential feed fed-batch 
[~, FBconst] = min(abs(tspan-16.3028));   %initiate constant feed fed-batch at 16.3028
[~, Tpulse1] = min(abs(tspan-20.5833));
[~, Tpulse2] = min(abs(tspan-22.3639));
[~, Tpulse3] = min(abs(tspan-24.4472));
[~, Tpulse4] = min(abs(tspan-29.5111));
[~, Tkla]    = min(abs(tspan-14.3417));
t_sample     = Offlinedata(:,1);
%% Initialization

y0 = [2 0.17 4.94 0.0129 98 98 0];   %Initial conditions
   %[V0, X0, S0, A0 DOT0 DOTm F0]
                                                        
Fe0 = 0.145*60/1000;             %L/h; initial value of exponential feed.

%Initial Parameter values (from literature)
      
        %[Kap     Ksa     Ko      Ks     Kia      Kis    pAmax   qAmax    qm     qSmax    Yas     Yoa     Yxa    Yem     Yos    Yxsof] 

Par = [0.5088	0.0128	0.0001	0.0381	1.2602	1.8383	0.2286	0.1148	0.0133	0.635	0.8938	0.5221	0.5794	0.5321	1.5722	0.229];

%% Inputs
Si = 300;                        % concentration of glucose feed 300 g/L
mufeed = 0.222;                  % set spec. growth rate during exp. feed
DOTstar = 99;                    % equil. DO concentration at operating pressure & temperature
Kla = 220;
tau = 35;                        % response time of DO probe

u = [Si mufeed DOTstar Kla tau]; %inputs vector


%%  Batch Phase
T = [];
Y = [];
flag1 = true; flag2 = false; flag3 = false; %specify which equations to solve using flags
tspan1 = tspan(1:FBinit,1);
options = odeset('NonNegative',1:7);
[t1, y1] = ode15s(@fn_e_coli,tspan1,y0,options,Par,u);
T = [T;t1];
Y = [Y;y1];
y0 = y1(end,:);

%% Exponential feed fed-batch
tspan21 = tspan(FBinit:Tkla,1);
flag1 = false; flag2 = true;
y0(7) = Fe0;
[t21, y21] = ode15s(@fn_e_coli,tspan21,y0,options,Par,u);

tspan22 = tspan(Tkla:FBconst,1);
u(4) = 355;                             %kla increased by inc. rpm
y0 = y21(end,:);  
[t22, y22] = ode15s(@fn_e_coli,tspan22,y0,options,Par,u);
t2 = vertcat(t21,t22(2:end));
y2 = vertcat(y21,y22(2:end,:));
%% Constant feed fed-batch, use half of the last value of exponential feed as
%   constant feed with intermittent pulses of glucose
y0 = y2(end,:);
y0(7) = 0.5*y0(7); 
tspan3 = tspan(FBconst:Tpulse1,1);
flag2 = false; flag3 = true;                           
[t3, y3] = ode15s(@fn_e_coli,tspan3,y0,options,Par,u);
y0 = y3(end,:);
y0(3) = 1.0; %give a pulse of 1g/L)
tspan4 = tspan(Tpulse1:Tpulse2,1);
[t4, y4] = ode15s(@fn_e_coli,tspan4,y0,options,Par,u);
y0 = y4(end,:);
y0(3) = 0.04; %give a pulse of 0.04g/L)
tspan5 = tspan(Tpulse2:Tpulse3,1);
[t5, y5] = ode15s(@fn_e_coli,tspan5,y0,options,Par,u);
y0 = y5(end,:);
y0(3) = 0.2; %give a pulse of 0.2g/L)
tspan6 = tspan(Tpulse3:Tpulse4,1);
[t6, y6] = ode15s(@fn_e_coli,tspan6,y0,options,Par,u);
y0 = y6(end,:);
y0(3) = 0.2; %give a pulse of 0.2g/L)
tspan7 = tspan(Tpulse4:end,1);
[t7, y7] = ode15s(@fn_e_coli,tspan7,y0,options,Par,u);
%% Collect all data points and extract corresponding offline data
T = [T;t2(2:end);t3(2:end);t4(2:end);t5(2:end);t6(2:end);t7(2:end)];
Y = [Y;y2(2:end,:);y3(2:end,:);y4(2:end,:);y5(2:end,:);y6(2:end,:);y7(2:end,:)];

%%  Plot profiles
FS = 12;
PW = 1.2;       %plotwidth: thickness of plotted line
BW = 2;         %thickness of boundaries

figure(2);
ax1 = subplot(2,2,1);
err = Offlinedata(:,3);
errorbar(t_sample,Offlinedata(:,2),err,'*r'); hold on
plot(T,Y(:,2),'b','LineWidth',PW);
plot([11.4469,11.4469],[0,25],'-.k');
plot([16.3028,16.3028],[0,25],'-.k');hold off
xlabel('Time(h)');ylabel('X (g/L)'); title('Biomass');
legend('Data','Model','Location','East'); legend('boxoff');
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

ax2 = subplot(2,2,2);
err = Offlinedata(:,5);
errorbar(t_sample,Offlinedata(:,4),err,'*r'); hold on
plot(T,Y(:,3),'b','LineWidth',PW);
plot([11.4469,11.4469],[0,5],'-.k');
plot([16.3028,16.3028],[0,5],'-.k');hold off
xlabel('Time(h)');ylabel('S (g/L)'); title('Glucose');
legend('Data','Model','Location','NorthEast');legend('boxoff');
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

ax3 = subplot(2,2,3);
err = Offlinedata(:,7);
errorbar(t_sample,Offlinedata(:,6),err,'*r'); hold on
plot(T,Y(:,4),'b','LineWidth',PW)
plot([11.4469,11.4469],[0,.35],'-.k');
plot([16.3028,16.3028],[0,.35],'-.k');hold off
xlabel('Time(h)');ylabel('A (g/L)'); title('Acetate');
legend('Data','Model');legend('boxoff');
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

ax4 = subplot(2,2,4);
plot(RawDOT(:,1),RawDOT(:,7),'r:','LineWidth',1.5); hold on
plot(T,Y(:,6),'b','LineWidth',PW)
plot([11.4469,11.4469],[0,100],'-.k');
plot([16.3028,16.3028],[0,100],'-.k');hold off
xlabel('Time(h)');ylabel('DOT (%)'); title('Dissolved Oxygen');
legend('Data','Model');legend('boxoff');
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

linkaxes([ax1,ax2,ax3,ax4],'x')
axis([0 35 -inf inf])

