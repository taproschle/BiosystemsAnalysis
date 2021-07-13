%% Comparison Plots
clear; clc;
% Data loading obtained from the simulation with the adjusted parameters

muset = 0.11;
DOset = 0.0024;
load 'simdata_xu.mat' dataXu; Xu = dataXu;
load 'simdata_dew.mat' dataDew; Dew = dataDew;
load 'simdata_an.mat' dataAn; An = dataAn;
Exp = table2array(readtable('data.csv'));

c1 = "#1B9E77"; c2 = "#D95F02"; c3 = "#7570B3"; c4 = "#E7298A";
c5 = "#66A61E"; c6 = "#E6AB02"; c7 = "#A6761D"; c8 = "#666666";

cAn  = [0.129 0.565 0.549];  lineAn  = '-';
cDew = [0.267 0.004 0.329];  lineDew = '-';
cXu  = [0.976 0.549 0.039];  lineXu  = '-';

%or2 = [0.988 0.631 0.031];
%% Biomass:
figure(1);
plot(Exp(:,1),Exp(:,2),'s','Color','k','linewidth',1); hold on;
title('Biomass Concentration','fontname','Calibri','fontsize',12)
plot(An(:,1),An(:,2),lineAn,'Color',cAn,'linewidth',1.5)
plot(Dew(:,1),Dew(:,2),lineDew,'Color',cDew,'linewidth',1.5)
plot(Xu(:,1),Xu(:,2),lineXu,'Color',cXu,'linewidth',1.5)
legend('Exp. Data','Anane','Dewasme','Xu','location','northwest',...
    'fontname','Calibri')
xlabel('Time [hr]')
ylabel('Biomass [gr/L]')
grid on
xlim([0 40.1])

%% Glucose
figure(2);
plot(Exp(:,1),Exp(:,3),'s','Color','k','linewidth',1); hold on;
title('Glucose Concentration','fontname','Calibri','fontsize',12)
plot(An(:,1),An(:,3),lineAn,'Color',cAn,'linewidth',1.5)
plot(Dew(:,1),Dew(:,3),lineDew,'Color',cDew,'linewidth',1.5)
plot(Xu(:,1),Xu(:,3),lineXu,'Color',cXu,'linewidth',1.5)
legend('Exp. Data','Anane','Dewasme','Xu','fontname','Calibri')
xlabel('Time [hr]')
ylabel('Glucose [gr/L]')
grid on
xlim([0 40.1])

%% Ethanol
figure(3);
plot(Exp(:,1),Exp(:,4),'s','Color','k','linewidth',1); hold on;
title('Ethanol Concentration','fontname','Calibri','fontsize',12)
plot(An(:,1),An(:,4),lineAn,'Color',cAn,'linewidth',1.5)
plot(Dew(:,1),Dew(:,4),lineDew,'Color',cDew,'linewidth',1.5)
plot(Xu(:,1),Xu(:,4),lineXu,'Color',cXu,'linewidth',1.5)
legend('Exp. Data','Anane','Dewasme','Xu','fontname','Calibri')
xlabel('Time [hr]')
ylabel('Ethanol [gr/L]')
grid on
xlim([0 40.1])

%% Reactor volume
figure(4);
plot(An(:,1),An(:,6),lineAn,'Color',cAn,'linewidth',1.5); hold on;
title('Reactor Volume','fontname','Calibri','fontsize',12)
plot(Dew(:,1),Dew(:,6),lineDew,'Color',cDew,'linewidth',1.5)
plot(Xu(:,1),Xu(:,6),'-','Color',cXu,'linewidth',1.5)
legend('Anane','Dewasme','Xu','location','northwest','fontname','Calibri')
xlabel('Time [hr]')
ylabel('Volume [L]')
grid on

%% Dissolved Oxygen
figure(5);
subplot(3,1,1)
plot(An(:,1),An(:,5),lineAn,'Color',cAn,'linewidth',1.5); hold on;
plot([0 Dew(end,1)],[DOset DOset],'--k','linewidth',1); hold off;
grid on
title('Dissolved Oxygen Concentration','fontname','Calibri','fontsize',12)
legend('Anane','DO_{set}','location','northwest','fontname','Calibri')
xlabel('Time [hr]')
ylabel('O_2 [gr/L]')
subplot(3,1,2)
plot(Dew(:,1),Dew(:,5),lineDew,'Color',cDew,'linewidth',1.5); hold on;
plot([0 Dew(end,1)],[DOset DOset],'--k','linewidth',1); hold off;
grid on
legend('Dewasme','DO_{set}','location','northwest','fontname','Calibri')
xlabel('Time [hr]')
ylabel('O_2 [gr/L]')
subplot(3,1,3)
plot(Xu(:,1),Xu(:,5),lineXu,'Color',cXu,'linewidth',1.5); hold on;
plot([0 Dew(end,1)],[DOset DOset],'--k','linewidth',1); hold off;
legend('Xu','DO_{set}','location','northwest','fontname','Calibri')
xlabel('Time [hr]')
ylabel('O_2 [gr/L]')
grid on

%% Growth Constants
figure(6);
%plot(An(:,1),An(:,7),lineAn,'Color',cAn,'linewidth',1.5); hold on;
plot(Dew(:,1),Dew(:,7),lineDew,'Color',cDew,'linewidth',1.5); hold on;
title('Growth constants','fontname','Calibri','fontsize',12)
%plot(Xu(:,1),Xu(:,7),lineXu,'Color',cXu,'linewidth',1.5)
plot([0 Dew(end,1)],[muset muset],'--k','linewidth',1); hold off;
%legend('Anane','Dewasme','Xu','\mu_{set}','fontname','Calibri')
legend('Dewasme','\mu_{set}','fontname','Calibri')
xlabel('Time [hr]')
ylabel('\mu [hr^{-1}]')
grid on

%% Critical growth and Critical substrate in dewasme
critical = 'y';
figure(7);
subplot(2,1,1)
plot(Dew(:,1),Dew(:,7),lineDew,'Color',cDew,'linewidth',1.5); hold on;
plot(Dew(:,1),Dew(:,8),lineDew,'Color',critical,'linewidth',1.5)
title('Critical growth in Dewasme Model','fontname','Calibri','fontsize',12)
legend('\mu','\mu_{crit}','location','best','fontname','Calibri')
xlabel('Time [hr]')
ylabel('\mu [hr^{-1}]')
grid on
subplot(2,1,2)
plot(Dew(:,1),Dew(:,3),lineDew,'Color',cDew,'linewidth',1.5); hold on;
plot(Dew(:,1),Dew(:,9),lineDew,'Color',critical,'linewidth',1.5)
title('Critical substrate in Dewasme Model','fontname','Calibri','fontsize',12)
legend('S','S_{crit}','location','best','fontname','Calibri')
xlabel('Time [hr]')
ylabel('S [gr/L]')
grid on

%% Inlet plots
figure(8);
subplot(2,2,1)
plot(An(:,1),An(:,7),lineAn,'Color',cAn,'linewidth',1.5); hold on;
plot(Dew(:,1),Dew(:,10),lineDew,'Color',cDew,'linewidth',1.5)
plot(Xu(:,1),Xu(:,7),lineXu,'Color',cXu,'linewidth',1.5); hold off
title('Inlet Flowrate','fontname','Calibri','fontsize',12)
legend('Anane','Dewasme','Xu','location','best','fontname','Calibri')
xlabel('Time [hr]')
ylabel('F [L/hr]')
grid on

subplot(2,2,2)
plot(An(:,1),An(:,8),lineAn,'Color',cAn,'linewidth',1); hold on;
plot(Dew(:,1),Dew(:,11),lineDew,'Color',cDew,'linewidth',1)
plot(Xu(:,1),Xu(:,8),lineXu,'Color',cXu,'linewidth',1); hold off
title('Agitation Rate','fontname','Calibri','fontsize',12)
legend('Anane','Dewasme','Xu','location','best','fontname','Calibri')
xlabel('Time [hr]')
ylabel('N [rpm]')
ylim([297 603])
grid on

subplot(2,2,3)
plot(An(:,1),An(:,9),lineAn,'Color',cAn,'linewidth',1); hold on;
plot(Dew(:,1),Dew(:,12),lineDew,'Color',cDew,'linewidth',1)
plot(Xu(:,1),Xu(:,9),lineXu,'Color',cXu,'linewidth',1); hold off
title('Aeration rate','fontname','Calibri','fontsize',12)
legend('Anane','Dewasme','Xu','location','best','fontname','Calibri')
xlabel('Time [hr]')
ylabel('G [L/min]')
ylim([0.195 0.605])
grid on

subplot(2,2,4)
plot(An(:,1),An(:,10),lineAn,'Color',cAn,'linewidth',1); hold on;
plot(Dew(:,1),Dew(:,13),lineDew,'Color',cDew,'linewidth',1)
plot(Xu(:,1),Xu(:,10),lineXu,'Color',cXu,'linewidth',1); hold off
title('O_2% of Inlet Airflow','fontname','Calibri','fontsize',12)
legend('Anane','Dewasme','Xu','location','best','fontname','Calibri')
xlabel('Time [hr]')
ylabel('yO_2 [L/min]')
ylim([0.2 1.01])
grid on




