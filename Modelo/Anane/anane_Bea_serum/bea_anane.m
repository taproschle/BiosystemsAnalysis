clear ; clc ; close all
% Este si es DAG
% En la linea 64 de "intento_anane.m" esta la ley de control
% Ahi modifica la weaita de acuerdo a lo que estan haciendo 

y0 = [6.8 0.4 0.1 0 0.75];   %Initial conditions
   %[V0, X0, S0, A0 DOT0]
                                                        
Fe0 = 0.145*60/1000;             %L/h; initial value of exponential feed.

%Initial Parameter values (from literature)
      
k   = [ 0.5088  % Kap
        0.0128  % Ksa
        0.0381  % Ks
        1.2602  % Kia
        1.8383  % Kis
        0.2286  % pAmax
        0.1148  % qAmax
        0.0133  % qm
        0.6350  % qSmax
        0.8938  % Yas
        0.5221  % Yoa
        0.5794  % Yxa
        0.5321  % Yem
        1.5722  % Yos
        0.229]; % Yxsof


% mufeed 0.18 optimo

Si = 550;                        % concentration of glucose feed 300 g/L
mufeed = 0.21;                    % set spec. growth rate during exp. feed
DOTstar = 0.85;                    % equil. DO concentration at operating pressure & temperature
Kla = 2200;
Ko  = 0.1;
tau = 35;                        % response time of DO probe

u = [Si mufeed DOTstar Kla Ko tau]; %inputs vector

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

