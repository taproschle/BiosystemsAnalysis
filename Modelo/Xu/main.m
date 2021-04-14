clear mex;
clear all;
close all;
global Vin Xin Factor_mu
%================================== DATA ==================================
load datos.txt;
texp  = datos(:,1)';
yexp = datos(:,2:end);
%====================== PROBLEM DEPENDENT DECLARATIONS=====================
Vin = 0.3; % [L]
Xin = yexp(1,1); % [g/L]
Factor_mu = 1.2;  
%=================== END OF PROBLEM DEPENDENT DECLARATIONS ================
%%
%========================= PROBLEM SPECIFICATIONS ===========================

problem.f='fun_costo'; %mfile containing the objective function

%parámetros = YS_ox_X  - YS_ox_O - YO_ox_X - YS_f_X - YS_f_IP - YIP_X -YIP_O - qS_max  - qIP_max - qO_max
% - K_S - K_IP - Ki_IP - K_O
problem.x_L = [ 0.5     400         1e-4   1e-3     1e-3      20      4         0.01     60     1e-3    1e-2       9     1e-5 ];
problem.x_U = [  1      800           0.1     0.5       0.6       80      6           5        100     0.3        5         20    0.01];
problem.x_0 = [0.68    483          0.07   0.46      0.28   39.58  5.6   0.2     84.15  0.23   0.35     15.5  0.005 ]';

opts.maxeval = 1.5e3; 
opts.maxtime = 1e4;
opts.local.n1  = 1000;
%opts.local.solver='fminsearch';
%opts.local.solver='solnp';
%opts.local.solver='fmincon';
%opts.local.solver='fminsearch';
opts.local.solver='IPOPT';
%opts.local.solver='misqp';

% opts.local.tol=3;
% %Try using this option and not using it (you'll see a great difference)
% opts.log_var=[1:8];              %Declare all variables as "logarithmic"
% opts.maxtime=1e6;
%========================= END OF PROBLEM SPECIFICATIONS ==================

%================================== OPTIMIZATION ==========================
Results=MEIGO(problem,opts,'ESS',texp,yexp);
 %Results=ssm_kernel(problem,opts,texp,yexp);
 
%%
%================================== Results ===============================
x  = Results.xbest;
X1 = [yexp(1,1) yexp(1,2) yexp(1,3) 0.3 7.3];   %condiciones inciciales
options = odeset('NonNegative',[1 2 3 4 5], 'RelTol', 1e-7, 'AbsTol', 1e-7);
%tspan = linspace(texp(1),texp(end),1e4);
[tout,yout] = ode15s(@model,texp,X1,options,x);

figure(1)
subplot(2,2,1); % Biomasa
plot(tout,yout(:,1),'-r',texp,yexp(:,1),'*b','LineWidth',1,'MarkerSize',3);
xlabel(' [h]'); ylabel('[g/L]'); grid on
title('\fontsize{14} Biomasa vs tiempo')

subplot(2,2,2); % Sustrato
plot(tout,yout(:,2),'-r',texp,yexp(:,2),'*b','LineWidth',1,'MarkerSize',3);
xlabel(' [h]'); ylabel('[g/L]'); grid on
title('\fontsize{14} Sustrato vs tiempo')

subplot(2,2,3); % Producto
plot(tout,yout(:,3),'-r',texp,yexp(:,3),'*b','LineWidth',1,'MarkerSize',3);
xlabel(' [h]'); ylabel('[g/L]'); grid on
title('\fontsize{14} Producto vs tiempo')

subplot(2,2,4); % Volumen
plot(tout,yout(:,4),'-r');
xlabel(' [h]'); ylabel('[L]'); grid on
title('\fontsize{14} Volumen vs tiempo')
