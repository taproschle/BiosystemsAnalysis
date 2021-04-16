%% Plots Anane
clear ; clc ; close all
c1 = [0.99607 0.850984 0.46274]';
c2 = [0.99215 0.552941 0.23529]';
c3 = [0.98823 0.305882 0.16470]';
c4 = [0.89019 0.101960 0.10980]';
reds = [c4 c3 c2 c1];
b1 = [0.61967 0.792157 0.88235]';
b2 = [0.41960 0.682353 0.83921]';
b3 = [0.25882 0.572549 0.77647]';
b4 = [0.12941 0.443137 0.70980]';
blues = [b1 b2 b3 b4];
g1 = [0.63137 0.850980 0.60784]';
g2 = [0.45902 0.768627 0.46274]';
g3 = [0.25490 0.670588 0.36470]';
g4 = [0.00000 0.352941 0.19607]';
greens = [g1 g2 g3 g4];
forms = ["-" , "--" , "-" , "--"];

% Este si es DAG
% En la linea 64 de "intento_anane.m" esta la ley de control
% Ahi modifica la weaita de acuerdo a lo que estan haciendo 

y0 = [6.8 0.4 0.1 0.0129 98 98];   %Initial conditions
   %[V0, X0, S0, A0 DOT0 DOTm F0]
                                                        
Fe0 = 0.145*60/1000;             %L/h; initial value of exponential feed.

%Initial Parameter values (from literature)
      
        %[Kap     Ksa     Ko      Ks     Kia      Kis    pAmax   qAmax    qm     qSmax    Yas     Yoa     Yxa    Yem     Yos    Yxsof] 

Par = [0.5088	0.0128	0.0001	0.0381	1.2602	1.8383	0.2286	0.1148	0.0133	0.635	0.8938	0.5221	0.5794	0.5321	1.5722	0.229];

% mufeed 0.18 optimo

Si = 350;                        % concentration of glucose feed 300 g/L
mufeed = 0.3;                    % set spec. growth rate during exp. feed
DOTstar = 99;                    % equil. DO concentration at operating pressure & temperature
Kla = 220;
tau = 35;                        % response time of DO probe

u1 = [Si 0.1 DOTstar Kla tau]';
u2 = [Si 0.15 DOTstar Kla tau]';
u3 = [Si 0.2 DOTstar Kla tau]';
u4 = [Si 0.25 DOTstar Kla tau]';
uT = [u1 u2 u3 u4];

time_span = [0 35];
options = odeset('NonNegative',1:6);

figure(1);
%subplot(3,1,1)
for i = 1:4
    [t1, Y] = ode15s(@intento_anane,time_span,y0,options,Par,uT(:,i)');
    plot(t1, Y(:,2) ,'Color',reds(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
legend('\mu_{set} = 0.1','\mu_{set} = 0.15','\mu_{set} = 0.2','\mu_{set} = 0.25','location','northwest','linewidth',0.75)
title('Effect of \mu_{set} on Biomass')
xlabel('Time [hr]')
ylabel('Biomass [g/l]')

figure(2);
subplot(2,1,1)
for i = 1:4
    [t1, Y] = ode15s(@intento_anane,time_span,y0,options,Par,uT(:,i)');
    plot(t1, Y(:,3) ,'Color',blues(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
legend('\mu_{set} = 0.1','\mu_{set} = 0.15','\mu_{set} = 0.2','\mu_{set} = 0.25','location','northwest','linewidth',0.75)
title('Effect of \mu_{set} on Substrate')
xlabel('Time [hr]')
ylabel('Substrate [g/l]')

subplot(2,1,2)
for i = 1:4
    [t1, Y] = ode15s(@intento_anane,time_span,y0,options,Par,uT(:,i)');
    plot(t1, Y(:,4) ,'Color',greens(:,i),'linewidth',1,'LineStyle',forms(i))
    hold on
end
legend('\mu_{set} = 0.1','\mu_{set} = 0.15','\mu_{set} = 0.2','\mu_{set} = 0.25','location','northwest','linewidth',0.75)
title('Effect of \mu_{set} on Acetate')
xlabel('Time [hr]')
ylabel('Acetate [g/l]')
