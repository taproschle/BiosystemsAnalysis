%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reg_analysis
% Performs al the regression analysis for a parameter set: significance,
% sensitivity and identifiability analysis, and BIC and AICc calculations.
% Returns all the different regression analysis outputs.
%
% Benjam�n J. S�nchez
% Last Update: 2014-07-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       [AICc,BIC,CI,CC,Mc,s,diff]
function [AICc,BIC,CI,CC] = reg_analysis_mod(k) %[AICc,BIC,CI,CC,Mc,sensib_prom,diff] = reg_analysis(k)

%Cargamos datos operacionales extras necesarios para este paso
Operational_Data = evalin('base','int_time');  int_time = Operational_Data(:,1);


%Load experimental data:
ydata    = evalin('base','ydata');
[~,xmod] = solve_ODE(k,int_time);  %IMPORTANTE: Si los estados observados no correponden a todos los estados simulados, se 
                                %debe reescribir Solve_ODE y corregir lsq_func para considerar solo los estados observados.
[~,yexp] = obj_var(xmod,ydata);                                         

%Jacobian and residual calculations using lsqcurvefit:
kL = log(exp(k).*0.9999999999);
kU = log(exp(k).*1.0000000001 + ones(1,length(k))*1e-15);
%save x k kL kU
[k2,~,res,~,~,~,Jac] = lsqcurvefit(@lsq_func,k,int_time,yexp,kL,kU,[]);

%AIC & BIC calculations:
SS    = sum(sum(res.^2));
m     = length(k2);
[N,~] = size(res);
AIC   = N*log(SS/N)+2*(m+1);
AICc  = AIC + 2*(m+1)*(m+2)/(N-m-2);
BIC   = SS/N + m*log(N);

%Confidence intervals:
[CI,CC] = intconfianza(Jac,res,k2,0.05);
CC      = CC./100;
disp('Insignificant Parameters:');
for i = 1:m
    if k2(i) == 0
        CC(i) = 0;
    elseif CC(i) >= 2
        disp(num2str(i));
    end
end

% %Identifiability analysis:
% x0 = evalin('base','x0');
% U  = evalin('base','U');
% Mc = identificaBSB(1:length(x0),k,x0,int_time,@model,U,0);
% set(0,'ShowHiddenHandles','on');delete(get(0,'Children'))

% %Sensitivity analysis:
% ksensibilidadBSB(1:length(x0),k,x0,texp,@model,3,0)
% sensib_prom = load('GraficoBarra.txt');
% set(0,'ShowHiddenHandles','on');delete(get(0,'Children'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%