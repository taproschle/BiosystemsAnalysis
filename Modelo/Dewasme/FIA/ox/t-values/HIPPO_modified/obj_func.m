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

function [J,g,R]=obj_func(k,~,ydata)

Operational_Data = evalin('base','int_time');
x0               = evalin('base','x0');

%Definicion del tiempo de integracion en base a datos operacionales de seguimiento de temperatura
int_time         = Operational_Data(:,1);

%Definimos dos tiempos de integracion: Previo a la adicion de FDA y posterior a ello.
t1        = int_time;

%Opciones de integracion
options = odeset('RelTol',1e-3,'AbsTol',1e-3);

%Integracion de la seccion previa a la adicion de FDA
[~,X]   = ode15s(@model,t1,x0,options,k);

ran = randi(3);
if ran == 1
    points = ".";
elseif ran == 2
    points = "..";
else
    points = "...";
end
disp('Iterating'+" "+points)

%Construction of residual matrix:
%Calculation of the objective function:

if ~ isequal(size(X(:,1:3)),size(ydata))
    J = 1e20;
    R = 1e10*ones(60,1);
else
    R =[(ydata(:,1)-X(:,1))./max(ydata(:,1)) (ydata(:,2)-X(:,2))./max(ydata(:,2)) ...
        (ydata(:,3)-X(:,3))./max(ydata(:,3))];
    J = sum(sum(R.^2));
    g = 0;
    R=reshape(R,numel(R),1);
end




    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%