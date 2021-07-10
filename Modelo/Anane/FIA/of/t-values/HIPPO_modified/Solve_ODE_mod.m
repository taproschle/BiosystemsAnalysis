function [T,Xf] = Solve_ODE_mod(k)    

Operational_Data = evalin('base','int_time');
x0               = evalin('base','x0');

%Definicion del tiempo de integracion en base a datos operacionales de seguimiento de temperatura
int_time         = Operational_Data(:,1);

%Error function definition for parameters estimation
t1        = int_time;

%Pre DAP addition section
[T,X]   = ode15s(@model,t1,x0,[],k);


%Resulting data consolidation
Xf        = X;
end