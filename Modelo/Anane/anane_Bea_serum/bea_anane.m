clear ; clc ; close all
% Este si es DAG
% En la linea 64 de "intento_anane.m" esta la ley de control
% Ahi modifica la weaita de acuerdo a lo que estan haciendo 

y0 = [6.8 0.4 0.1 0 0.75];   %Initial conditions
   %[V0, X0, S0, A0 DOT0]
                                                        
Fe0 = 0.145*60/1000;             %L/h; initial value of exponential feed.

%Initial Parameter values (from literature)
      
        %[Kap     Ksa     Ko      Ks     Kia      Kis    pAmax   qAmax    qm     qSmax    Yas     Yoa     Yxa    Yem     Yos    Yxsof] 

Par = [0.5088	0.0128	0.0001	0.0381	1.2602	1.8383	0.2286	0.1148	0.0133	0.635	0.8938	0.5221	0.5794	0.5321	1.5722	0.229];

% mufeed 0.18 optimo

Si = 550;                        % concentration of glucose feed 300 g/L
mufeed = 0.21;                    % set spec. growth rate during exp. feed
DOTstar = 0.85;                    % equil. DO concentration at operating pressure & temperature
Kla = 2200;
tau = 35;                        % response time of DO probe

u = [Si mufeed DOTstar Kla tau]; %inputs vector

time_span = [0 35];

options = odeset('NonNegative',1:5);
[t1, y1] = ode15s(@anane_bea,time_span,y0,options,Par,u);

figure(1) %V X S A
subplot(2,3,1)
plot(t1,y1(:,1))
legend('V')
subplot(2,3,2)
plot(t1,y1(:,2))
legend('X')
subplot(2,3,3)
plot(t1,y1(:,3))
legend('S')
subplot(2,3,4)
plot(t1,y1(:,4))
legend('A')
subplot(2,3,5)
plot(t1,y1(:,5))
legend('O2')


% x0=10;
% y0=10;
% width=900;
% height=600;
% set(fig,'position',[x0,y0,width,height])

%% 
