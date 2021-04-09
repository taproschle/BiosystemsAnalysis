fig = figure(1);
FS = 12; FontWeight = 'bold'; % Fontsize 15 and bold for the plots
off_idx        = [1   6   7   18    19  30  31  33];  %indices of offline data X = [1 9], S = [10 23]...

ax1 = subplot(3,2,1);
plot(T,Y(:,1),'b'); hold on
plot(Offlinedata1(off_idx(1):off_idx(2),1),Offlinedata1(off_idx(1):off_idx(2),2),'r*'); hold off
ylabel('X (g/L)'); xlabel('Fermentation Time (h)')
title('Biomass')
set(gca,'LineWidth',2,'FontSize',FS,'FontWeight',FontWeight) 
legend('Model','Data','Location','NorthWest');legend('boxoff')


ax2 = subplot(3,2,2);
plot(T,Y(:,2),'b'); hold on
plot(Offlinedata1(off_idx(3):off_idx(4),1),Offlinedata1(off_idx(3):off_idx(4),2),'r*'); hold off
ylabel('S (g/L)'); xlabel('Fermentation Time (h)')
title('Glucose')
set(gca,'LineWidth',2,'FontSize',FS,'FontWeight',FontWeight) 
legend('Model','Data','Location','NorthEast');legend('boxoff')


ax3 = subplot(3,2,3);
plot(T,Y(:,3),'b'); hold on
plot(Offlinedata1(off_idx(5):off_idx(6),1),Offlinedata1(off_idx(5):off_idx(6),2),'r*'); hold off
ylabel('A (g/L)','Fontsize',8); xlabel('Fermentation Time (h)','Fontsize',8)
title('Acetate')
set(gca,'LineWidth',2,'FontSize',FS,'FontWeight',FontWeight) 
legend('Model','Data','Location','NorthEast');legend('boxoff')



ax4 = subplot(3,2,[4 6]);
plot(T, Y(:,4),'b'); hold on
plot(Fermdata1(FBinit:6879,1),Fermdata1(FBinit:6879,2),'r-.'); hold off
ylabel('DOT (%)'); xlabel('Fermentation Time (h)')
title('Dissolved Oxygen')
set(gca,'LineWidth',2,'FontSize',FS,'FontWeight',FontWeight) 
legend('Model','Data','Location','NorthEast');legend('boxoff')
 

ax5 = subplot(3,2,5);
plot(T,Y(:,5),'b'); hold on
plot(Offlinedata1(off_idx(7):off_idx(8),1),Offlinedata1(off_idx(7):off_idx(8),2),'r*'); hold off
ylabel('P (g/g)'); xlabel('Fermentation Time (h)')
title('Product')
set(gca,'LineWidth',2,'FontSize',FS,'FontWeight',FontWeight) 
legend('Model','Data','Location','NorthWest');legend('boxoff')

% linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
% axis([14 20 -inf inf])
% fig.PaperType = 'a4'; fig.PaperOrientation = 'landscape';%fig.WindowStyle = 'docked';


