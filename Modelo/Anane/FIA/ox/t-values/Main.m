clear ; clc ; close all

global int_time x0

%Load data and integration options
[kL,k0,kU,opts_SSm,texp,...
    ydata,x0,solver_ODE,opts_ODE,int_time] = load_problem;

id = 'MATLAB:ode15s:IntegrationTolNotMet';
warning('off',id)

%Select parameters to be fixed (manually)

fixed_params = [12 8 6 4 9 1 5]; %Indexes of parameters to be fixed 

kfixed = nan(12,1);
for k = fixed_params
    kfixed(k) = k0(k);
end

%Calculate optimal fitting parameters, CC and CI
[CI,CC,k_SSm] = iteration_mod(kfixed);

%Calculate t-values for this model structure
t_values = 4*k_SSm'./CI;