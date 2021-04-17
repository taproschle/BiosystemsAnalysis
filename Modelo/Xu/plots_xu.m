%% Plots Xu
clear
primerasfiguras
colors
global mu_set;
mu_set = 0.4;

tspan = [0 25];
% State variables plot
figure(3);
fig = figure(3);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
[T,C] = ode15s(@fedbatch, tspan, x0);
yyaxis left
plot(T,C(:,3),'-r','linewidth',1)
grid on
ylabel('X [g/L]')
hold on
yyaxis right
plot(T,C(:,1),'-b','linewidth',1)
plot(T,C(:,2),'-g','linewidth',1)
hold off
legend('X','S','P','location','east')
title('State variables for Xu Model')
xlabel('Time [hr]')
ylabel('S,P [g/L]')


%%

mu = [0.2 0.3 0.4 0.5];
figure(1);
for i = 1:4 
    mu_set = mu(i);
    [T,C] = ode15s(@fedbatch, tspan, x0);
    plot(T,C(:,3),'Color',reds(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
hold off
legend('\mu_{set} = 0.2','\mu_{set} = 0.3','\mu_{set} = 0.4','\mu_{set} = 0.5','location','northwest','linewidth',0.75)
title('Effect of \mu_{set} on Biomass (Xu Model)')
xlabel('Time [hr]')
ylabel('Biomass [g/l]')
grid on

figure(2);
subplot(2,1,1)
for i = 1:4 
    mu_set = mu(i);
    [T,C] = ode15s(@fedbatch, tspan, x0);
    plot(T,C(:,1),'Color',blues(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
hold off
legend('\mu_{set} = 0.2','\mu_{set} = 0.3','\mu_{set} = 0.4','\mu_{set} = 0.5','location','northwest','linewidth',0.75)
title('Effect of \mu_{set} on Substrate (Xu Model)')
xlabel('Time [hr]')
ylabel('Substrate [g/l]')
grid on

subplot(2,1,2)
for i = 1:4 
    mu_set = mu(i);
    [T,C] = ode15s(@fedbatch, tspan, x0);
    plot(T,C(:,2),'Color',greens(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
hold off
legend('\mu_{set} = 0.2','\mu_{set} = 0.3','\mu_{set} = 0.4','\mu_{set} = 0.5','location','northwest','linewidth',0.75)
title('Effect of \mu_{set} on Acetate (Xu Model)')
xlabel('Time [hr]')
ylabel('Acetate [g/l]')
grid on
