BW  = 1.5;
FS  = 16;


figure;
l = subplot(2,3,1);
plot(X_residuals,'ro');hold on
plot([0 20],[0 0],'k');hold off
xlabel('Sample Number'); ylabel('Y^{meas}-Y^{sim}');
title(l,'Residuals of Biomass')
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

l = subplot(2,3,2);
plot(S_residuals,'ro');hold on
plot([0 20],[0 0],'k');hold off
xlabel('Sample Number'); ylabel('Y^{meas}-Y^{sim}');
title(l,'Residuals of Substrate')
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

l = subplot(2,3,3);
plot(A_residuals,'ro');hold on
plot([0 20],[0 0],'k');hold off
xlabel('Sample Number'); ylabel('Y^{meas}-Y^{sim}');
title(l,'Residuals of Acetate')
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

l = subplot(2,3,4);
plot(DOT_residuals,'ro');hold on
plot([0 20],[0 0],'k');hold off
xlabel('Sample Number'); ylabel('Y^{meas}-Y^{sim}');
title(l,'Residuals of DOT')
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

l = subplot(2,3,5);
plot(All_res,'ro');hold on
plot([0 80],[0 0],'k');hold off
xlabel('Sample Number'); ylabel('Y^{meas}-Y^{sim}');
title(l,'All residuals pooled together')
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 

l = subplot(2,3,6);
plot(Tot_res,'ro');hold on
plot([0 20],[0 0],'k');hold off
xlabel('Sample Number'); ylabel('Y^{meas}-Y^{sim}');
title(l,'Summed residuals at sampling times')
set(gca,'LineWidth',BW,'FontSize',FS,'FontWeight','bold') 
