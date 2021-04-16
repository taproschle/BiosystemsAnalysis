clear ; clc ; close all
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

u = [Si mufeed DOTstar Kla tau]; %inputs vector

time_span = [0 35];

options = odeset('NonNegative',1:6);
[t1, y1] = ode15s(@intento_anane,time_span,y0,options,Par,u);

% Plot
fig = figure(1);
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
yyaxis left
plot(t1, y1(:, 2),'-r', "LineWidth", 1.2)
ylabel('Biomass [g/l]')
yyaxis right
ylabel('S and A [g/l]')
hold on
plot(t1, y1(:, 3),'-b' , "LineWidth", 1.2)
plot(t1, y1(:, 4),'-g' ,"LineWidth", 1.2)
grid on
legend("Biomass", "Sustrate", "Acetate", "LineWidth", 0.7, "location", "northwest")
title("Evolution of the system for E.coli")
xlabel("Time [h]")
hold off 
%% 
plot(t1, y1(:,5))