%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obj_func
% Integrates the dynamic model and calculates afterwards the cuadratic 
% difference between the model predictions and the experimental data,
% normalized by the maximum value of the experimental data. Returns the
% cuadratic difference (objective function). To be used in the parameter
% estimation with SSm.
%
% Benjam�n J. S�nchez
% Last Update: 2013-06-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J,g,R]=obj_func(k,texp,ydata)

%Resolution of ODEs:
[~,xmod] = solve_ODE(k,texp);

%Construction of residual matrix:
[ymod,yexp] = obj_var(xmod,ydata);
R = ymod - yexp;

%Calculation of the objective function:
J = sum(sum(R.^2));
disp(['J = ' num2str(J)]);

R=reshape(R,numel(R),1);
g=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%