fig = figure;
FS = 12; FontWeight = 'bold'; % Fontsize 15 and bold for the plots
% off_idx        = [1  33  34  66  67  99  100  123];  %indices of offline data X = [1 9], S = [10 23]...
Plot_data = data.Plot_data;
ax1 = subplot(3,2,1);
plot(T,Y(:,1),'b'); hold on
% plot(data.Offlinedata(off_idx(1):off_idx(2),1),0.37*(data.Offlinedata(off_idx(1):off_idx(2),2)),'r*'); hold off
plot(Plot_data(1).biomass(:,1),0.37*Plot_data(1).biomass(:,2),'r*')
plot(Plot_data(2).biomass(:,1),0.37*Plot_data(2).biomass(:,2),'k*')
plot(Plot_data(3).biomass(:,1),0.37*Plot_data(3).biomass(:,2),'c*'); hold off

ylabel('X (g/L)'); %xlabel('Fermentation Time (h)')
title('Biomass')
set(gca,'LineWidth',2,'FontSize',FS,'FontWeight',FontWeight, 'xlabel',[]) 
legend('Model','Data','Location','NorthWest');legend('boxoff')


ax2 = subplot(3,2,2);
plot(T,Y(:,2),'b'); hold on
% plot(data.Offlinedata(off_idx(3):off_idx(4),1),data.Offlinedata(off_idx(3):off_idx(4),2),'r*'); hold off
plot(Plot_data(1).Glu(:,1),Plot_data(1).Glu(:,2),'r*')
plot(Plot_data(2).Glu(:,1),Plot_data(2).Glu(:,2),'k*')
plot(Plot_data(3).Glu(:,1),Plot_data(3).Glu(:,2),'c*'); hold off

ylabel('S (g/L)');% xlabel('Fermentation Time (h)')
title('Glucose')
set(gca,'LineWidth',2,'FontSize',FS,'FontWeight',FontWeight) 
legend('Model','Data','Location','NorthEast');legend('boxoff')


ax3 = subplot(3,2,3);
plot(T,Y(:,3),'b'); hold on
% plot(data.Offlinedata(off_idx(5):off_idx(6),1),data.Offlinedata(off_idx(5):off_idx(6),2),'r*'); hold off
plot(Plot_data(1).Ace(:,1),Plot_data(1).Ace(:,2),'r*')
plot(Plot_data(2).Ace(:,1),Plot_data(2).Ace(:,2),'k*')
plot(Plot_data(3).Ace(:,1),Plot_data(3).Ace(:,2),'c*'); hold off

ylabel('A (g/L)','Fontsize',8); %xlabel('Fermentation Time (h)','Fontsize',8)
title('Acetate')
set(gca,'LineWidth',2,'FontSize',FS,'FontWeight',FontWeight) 
legend('Model','Data','Location','NorthEast');legend('boxoff')



ax4 = subplot(3,2,[4 6]);
plot(T, Y(:,4),'b'); hold on
plot(tspan(data.kla_change(6):end),DOT(data.kla_change(6):end),'r-.'); hold off
ylabel('DOT (%)'); xlabel('Fermentation Time (h)')
title('Dissolved Oxygen')
set(gca,'LineWidth',2,'FontSize',FS,'FontWeight',FontWeight) 
legend('Model','Data','Location','NorthEast');legend('boxoff')
 

ax5 = subplot(3,2,5);
plot(T,1000.*Y(:,5),'b'); hold on
% plot(data.Offlinedata(off_idx(7):off_idx(8),1),data.Offlinedata(off_idx(7):off_idx(8),2),'r*'); hold off
plot(Plot_data(1).Pro(:,1),Plot_data(1).Pro(:,2),'r*')
plot(Plot_data(2).Pro(:,1),Plot_data(2).Pro(:,2),'k*')
plot(Plot_data(3).Pro(:,1),Plot_data(3).Pro(:,2),'c*'); hold off

ylabel('P (mg/g)'); xlabel('Fermentation Time (h)')
title('Product')
set(gca,'LineWidth',2,'FontSize',FS,'FontWeight',FontWeight) 
legend('Model','Data','Location','NorthWest');legend('boxoff')

% linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
% axis([14 20 -inf inf])
% fig.PaperType = 'a4'; fig.PaperOrientation = 'landscape';%fig.WindowStyle = 'docked';


