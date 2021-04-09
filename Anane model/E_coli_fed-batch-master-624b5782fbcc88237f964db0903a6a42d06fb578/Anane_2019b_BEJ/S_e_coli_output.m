function [T_out, Y_output,Jac]   = S_e_coli_output(tspan,Par,y0,Sens_calc)

%Emmanuel Anane
%28/07/2016

%Simulation of fed-batch culture of E. coli. The script consists of a 
%simple batch phase, followed by an exponential fed-batch (with no product
%formation) and then a constant feed fed-batch phase.
%Load laboratory data



%% Define time spans
[~, FBinit] = min(abs(tspan-11.4469));   %batch phase ends at 11.4469h, initiate exponential feed fed-batch 
[~, FBconst] = min(abs(tspan-16.3028));  %initiate constant feed fed-batch at 16.3028
[~, Tpulse1] = min(abs(tspan-20.5833));
[~, Tpulse2] = min(abs(tspan-22.3639));
[~, Tpulse3] = min(abs(tspan-24.4472));
[~, Tpulse4] = min(abs(tspan-29.5111));
[~, Tkla] =  min(abs(tspan-14.3417));
%% Initialization

% y0 = [2 0.17 4.94 0.0129 98 98 0];   %Initial conditions
   %[V0, X0, S0, A0 DOT0 DOTm F0]
Number_state_var = 7;                                                        
Fe0 = 0.145*60/1000;             %L/h; initial value of exponential feed.

%Initial Parameter values (from literature)
      
        %[Kap     Ksa    Ks     Kia      Kis    pAmax   qAmax    qm     qSmax    Yas     Yoa     Yxa    Yem     Yos    Yxsof] 

data.Par = Par;
% Par =[0.505208595,0.013412432,0.0001,0.036997403,1.239949991,2.123050851,0.226768867000000,0.114762576000000,0.0129435980000000,0.635593525000000,0.909716720000000,0.543977350000000,0.571839309000000,0.533250672000000,1.56195758200000,0.226813106000000];

% Par = [0.5088	0.0128	0.0001	0.0381	20.602	1.8383	0.2286	0.1148	0.0133	0.635	0.8938	0.5221	0.5794	0.5321	1.5722	0.229];
%% Inputs
Si = 300;                        % concentration of glucose feed 200 g/L
mufeed = 0.222;                  % set spec. growth rate during exp. feed
DOTstar = 99;                    % equil. DO concentration at operating pressure & temperature
Kla = 220;
tau = 35;
data.u = [Si mufeed DOTstar Kla tau]; %imputs vector

%%  Batch Phase
T = [];
Y = [];
data.flag = true; %specify which equations to solve using flags
tspan1 = tspan(1:FBinit,1);
options = odeset('NonNegative',1:7);
[t1, y1] = ode15s(@fn_e_coli,tspan1,y0,options,data);
T = [T;t1];
Y = [Y;y1];
y0 = y1(end,:);

%% Exponential feed fed-batch
tspan21 = tspan(FBinit:Tkla,1);
data.flag = false; 
y0(7) = Fe0;
[t21, y21] = ode15s(@fn_e_coli,tspan21,y0,options,data);

tspan22 = tspan(Tkla:FBconst,1);
data.u(4) = 355;                             %kla increased by inc. rpm
y0 = y21(end,:);  
[t22, y22] = ode15s(@fn_e_coli,tspan22,y0,options,data);
t2 = vertcat(t21,t22(2:end));
y2 = vertcat(y21,y22(2:end,:));
%% Constant feed fed-batch, use half of the last value of exponential feed as
%   constant feed. Give intermittent pulses of glucose
data.flag = true;
y0 = y2(end,:);
y0(7) = 0.5*y0(7); 
tspan3 = tspan(FBconst:Tpulse1,1);
[t3, y3] = ode15s(@fn_e_coli,tspan3,y0,options,data);
y0 = y3(end,:);
y0(3) = 1.0; %give a pulse of 1g/L)
tspan4 = tspan(Tpulse1:Tpulse2,1);
[t4, y4] = ode15s(@fn_e_coli,tspan4,y0,options,data);
y0 = y4(end,:);
y0(3) = 0.04; %give a pulse of 0.02g/L)
tspan5 = tspan(Tpulse2:Tpulse3,1);
[t5, y5] = ode15s(@fn_e_coli,tspan5,y0,options,data);
y0 = y5(end,:);
y0(3) = 0.2; %give a pulse of 0.1g/L)
tspan6 = tspan(Tpulse3:Tpulse4,1);
[t6, y6] = ode15s(@fn_e_coli,tspan6,y0,options,data);
y0 = y6(end,:);
y0(3) = 0.2; %give a pulse of 0.08g/L)
tspan7 = tspan(Tpulse4:end,1);
[t7, y7] = ode15s(@fn_e_coli,tspan7,y0,options,data);
%% Collect all data points and extract corresponding offline data
T_out = [T;t2(2:end);t3(2:end);t4(2:end);t5(2:end);t6(2:end);t7(2:end)];
Y_output = [Y;y2(2:end,:);y3(2:end,:);y4(2:end,:);y5(2:end,:);y6(2:end,:);y7(2:end,:)];


%% Extract absolute sensitivities from the simulated data 

if strcmp(Sens_calc,'yes')
        J                   = Y_output(:,(Number_state_var+1):end);  %Extract Jacobians and reshape to [1:1:t(end)*length(y0) x length(Par)] matrix

        [q, ~]              = size(J);
        Jac                 = zeros(q,length(Par));

            for k2          = 1:q
                R           = reshape(J(k2,:),[Number_state_var length(Par)]);
                S           = sum(R);
                Jac(k2,:)   = S;
            end
else
    Jac = 0;
end

end