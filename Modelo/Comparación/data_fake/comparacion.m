clear ; clc ; close all

load Andata
load Dewdata
load Xudata
load data.csv

texp    = data(:,1);
yexp    = data(:,2:6);

tAn     = Andata(:,1);
cAn     = Andata(:,2:6);

tDew    = Dewdata(:,1);
cDew    = Dewdata(:,2:6);

tXu     = Xudata(:,1);
cXu     = Xudata(:,2:6);

c1 = "#1B9E77"; c2 = "#D95F02"; c3 = "#7570B3"; c4 = "#E7298A";
c5 = "#66A61E"; c6 = "#E6AB02"; c7 = "#A6761D"; c8 = "#666666";

tiledlayout(2,3)

nexttile
plot(texp,yexp(:,1),'s','Color',c1,'LineWidth',1)
grid on
hold on
plot(tAn,cAn(:,1),'-.','Color',c1,'LineWidth',1.5)
plot(tDew,cDew(:,1),'-','Color',c1,'LineWidth',1.5)
plot(tXu,cXu(:,1),'--','Color',c1,'LineWidth',1.5)
ylabel('Biomass (g/L)')
xlabel('Time (h)')
legend('Data','Anane','Dewasme','Xu','Location','Best')

nexttile
plot(texp,yexp(:,2),'s','Color',c2,'LineWidth',1)
grid on
hold on
plot(tAn,cAn(:,2),'-.','Color',c2,'LineWidth',1.5)
plot(tDew,cDew(:,2),'-','Color',c2,'LineWidth',1.5)
plot(tXu,cXu(:,2),'--','Color',c2,'LineWidth',1.5)
ylabel('Glucose (g/L)')
xlabel('Time (h)')
legend('Data','Anane','Dewasme','Xu','Location','Best')

nexttile
plot(texp,yexp(:,3),'s','Color',c3,'LineWidth',1)
grid on
hold on
plot(tAn,cAn(:,3),'-.','Color',c3,'LineWidth',1.5)
plot(tDew,cDew(:,3),'-','Color',c3,'LineWidth',1.5)
plot(tXu,cXu(:,3),'--','Color',c3,'LineWidth',1.5)
ylabel('Ethanol (g/L)')
xlabel('Time (h)')
legend('Data','Anane','Dewasme','Xu','Location','Best')

nexttile
plot(texp,yexp(:,4),'s','Color',c4,'LineWidth',1)
grid on
hold on
plot(tAn,cAn(:,4),'-.','Color',c4,'LineWidth',1.5)
plot(tDew,cDew(:,4),'-','Color',c4,'LineWidth',1.5)
plot(tXu,cXu(:,4),'--','Color',c4,'LineWidth',1.5)
ylabel('Dissolved Oxygen (g/L)')
xlabel('Time (h)')
legend('Data','Anane','Dewasme','Xu','Location','Best')

nexttile
plot(texp,yexp(:,5),'s','Color',c5,'LineWidth',1)
grid on
hold on
plot(tAn,cAn(:,5),'-.','Color',c5,'LineWidth',1.5)
plot(tDew,cDew(:,5),'-','Color',c5,'LineWidth',1.5)
plot(tXu,cXu(:,5),'--','Color',c5,'LineWidth',1.5)
ylabel('Volume (L)')
xlabel('Time (h)')
legend('Data','Anane','Dewasme','Xu','Location','Best')






