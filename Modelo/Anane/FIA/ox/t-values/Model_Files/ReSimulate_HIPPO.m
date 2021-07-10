function [T,Xf] = ReSimulate_HIPPO(k)    
global int_time x0

%Error function definition for parameters estimation
t1        = int_time;

%Pre DAP addition section
[T,Xf]   = ode23s(@model,t1,x0,[],k);

end