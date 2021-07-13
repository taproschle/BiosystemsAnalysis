%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obj_func
% Integrates the dynamic model and calculates afterwards the cuadratic 
% difference between the model predictions and the experimental data,
% normalized by the maximum value of the experimental data. Returns the
% cuadratic difference (objective function). To be used in the parameter
% estimation with SSm.
%
% Benjamín J. Sánchez
% Last Update: 2013-06-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J,g,R]=obj_func(k,texp,ydata)

%Resolution of ODEs:
[~,xmod] = solve_ODE(k,texp);

if ~ isequal(size(xmod(:,1:3)),size(ydata))
    J = 1e20;
    R = 1e10*ones(14,1);
    g = NaN;
%     disp(['J =' num2str(J)]);
else
    R =[(ydata(:,1)-xmod(:,1))./max(ydata(:,1)) (ydata(:,2)-xmod(:,2))./max(ydata(:,2)) ...
        (ydata(:,3)-xmod(:,3))./max(ydata(:,3))];
    J = sum(sum(R.^2));
    g = 0;
    R=reshape(R,numel(R),1);
%  disp(['J =' num2str(J)]);
end

%Construction of residual matrix:
% [ymod,yexp] = obj_var(xmod,ydata);
% R = ymod - yexp;

%Calculation of the objective function:
% J = sum(sum(R.^2));
% disp(['J = ' num2str(J)]);
% 
% R=reshape(R,numel(R),1);
% g=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%