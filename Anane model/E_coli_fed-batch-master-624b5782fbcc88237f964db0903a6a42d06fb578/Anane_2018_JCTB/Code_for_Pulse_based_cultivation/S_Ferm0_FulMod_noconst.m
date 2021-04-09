%{ 

%Emmanuel Anane, 
%21/09/2017


Description: Simulation of fed-batch culture of E. coli and comparing model
outputs to experimental data of the pulse-based cultivation after model 
fitting. The script implements a pulse feeding profile after an initial 
batch phase as follows: 
(i) FBexp--Exponential fed-batch where the feed is switched on for 1 min,
and'off' for 9 minutes. During the 1 min, the mass of glucose delivered is
equal to what would have been fed exponentially during the next 9 minutes 
when the feed is off. This creates glucose concentration gradients. During
the feeding phase, the rpm is reduced significantly (lower kla) to subject
the culture to high glucose-low oxygen zones. When the feed is off, rpm
(kla) is increased again to create low glucose-high oxygen zones.
(ii) FBconst-- same feeding scheme as in the exponential feeding phase
except that the feed rate is constant. Protein production was induced in
this phase.
%}
%% Load experimental data
load Fermdata1;load Offlinedata1;load Pulse 

tspan               = Fermdata1(1:6879,1);
[~, FBinit]         = min(abs(tspan-14.17172222));   %batch phase ends at 14.1717h, initiate exponential feed fed-batch 
t_sample            = Offlinedata1(:,1);             % extract sampling times from offline data
%% Inputs (Experimental Conditions)

Si                  = 200;
DOTstar             = 99.8;
Kla                 = 350;
mu_set              = 0.2494;
tau                 = 10;
gamma               = 0.8;                             %strength of induction
F                   = 0;
V                   = 2;
Fe0                 = 0.235*60/1000;         %L/h; initial value of exp. feed.                                       %the equivalent pulse feed is administered.

data.u              = [Si DOTstar Kla tau gamma F V mu_set];
data.Kla_vector     = [350 250 650];
   
%% Initialization

deltsp              = 1/60;
tlog                = 10/3600;              %computer logs data every 10s
y0          = [2.68 0   0  77.2 0];

%------------- Parameter values from PE -----------------------
              %Kap      Ksa      Ks      Kia     Kis     pAmax       qAmax    qm       qSmax      Yas     Yoa      Yxa         Yem      Yxsof     Yos     Ypx 
data.Par       = [0.359999996146146,0.0200000020553886,0.0199999997430764,14.2689997618575,9.99999976876878,0.326200001896012,0.250000003853789,0.0400000012846179,0.519999998201535,0.970000000513598,1.19999998972306,0.500000010276815,0.481000003057267,0.849999980730732,1.07999999280614,0.531999999950806];


T = [];
Y = [];

%% Exponential Pulse fed-batch

data.flag   = 1; 

cycles      = 17;      %there were 17 pulses during exponential feeding

for k1      = 1:cycles
    T_start = Pulse(1,k1);               %Time to start the pulse
    deltse = Pulse(2,k1);                %Time span to calculate the exp. feed.   
    mu_set = data.u(8);
    tse = [T_start T_start+deltse];      %Integrate Fe in this interval(10min).
    tsp = tse(1):tlog:tse(1)+deltsp;     %Feed equival.Fp in this interval(1min).
    fe = Fe0*exp(mu_set*deltse);         %Analytical solution of F. 
    Intfe = (1./mu_set).*Fe0.*exp(mu_set*deltse)-...
    (1./mu_set)*Fe0*exp(mu_set*0);       %integrate to find area under F.
    fp = Intfe./deltsp;                  %define fp as area/time (L/h)
                                         %F = Par.Fp in the function to 
                                         %include Fp in the ODE parameters.

    if k1 == 1
        data.u(3)    = data.Kla_vector(1);             % there wasn't a Kla change during the first pulse.
    else
        data.u(3)    = data.Kla_vector(2);             % set Kla during feeding
    end

    data.u(6)    = fp;                                 % set pulse feed rate;

    V            = data.u(7)+(fp*deltsp);              % re-calculate volume.
    data.u(7)    = V;

    options      = odeset('NonNegative',1:5);
   [t2, y2]      = ode15s(@fn_Ecoli,tsp,y0,options,data);


    y0           = y2(end,:);
    
    if k1 <= 2
        data.u(3)    = data.Kla_vector(1);              % there wasn't a Kla change during the first pulse, 2nd pulse went up to 800rpm, as in batch phase.
    else
        data.u(3)    = data.Kla_vector(3);              % set Kla during resting phase   
    end
 
   data.u(6)   = 0;                                     % feed is switched off for 9 min
   tswof       = tsp(end):tlog:tse(end);                %tswof = time span without feed(9 min).
   [t3, y3]     = ode15s(@fn_Ecoli,tswof,y0,options,data);

    y0          = y3(end,:);                            %set new parameters for looping.
    Fe0         = fe;
    
    if k1   == 1
    T           = [T;t2(1:end);t3(2:end)];              %concatenate t's and y's.
    Y           = [Y;y2(1:end,:);y3(2:end,:)]; 
    else
    T           = [T;t2(2:end);t3(2:end)];              %concatenate t's and y's.
    Y           = [Y;y2(2:end,:);y3(2:end,:)]; 
    end

end
%% Constant pulse fed-batch
data.flag   = 2; 
cycles      = 13;                         %there were 13 pulses during the constant feed phase, with protein production.

for k2 = 1:cycles
    
    w       = k1+k2;                   
    T_start = Pulse(1,w);
    deltse  = Pulse(2,w);
    tse     = [T_start T_start+deltse];     
    tsp     = tse(1):tlog:tse(1)+deltsp;     
    
    data.u(8) = 0.1;
    mu_set = data.u(8);
    tse = [T_start T_start+deltse];      %Integrate Fe in this interval(10min).
    fe = Fe0*exp(mu_set*deltse);         %Analytical solution of F. 
    Intfe = (1./mu_set).*Fe0.*exp(mu_set*deltse)-...
    (1./mu_set)*Fe0*exp(mu_set*0);       %integrate to find area under F.
    fp = Intfe./deltsp;                  %define fp as area/time (L/h)
                                         %F = Par.Fp in the function to 
                                         %include Fp in the ODE parameters.

    
    data.u(6)   = fp;                               % set constant pulse feed rate;

    data.u(3)   = data.Kla_vector(2);                              % set Kla during feeding

    V           = data.u(7)+(data.u(6)*deltsp);     % re-calculate volume.
    data.u(7)   = V;

    [t4, y4]    = ode15s(@fn_Ecoli,tsp,y0,options,data);

    y0          = y4(end,:);
    tswof       = tsp(end):tlog:tse(end); 
    data.u(3)   = data.Kla_vector(3);
    data.u(6)   = 0;                               % feed is switched off for 9 min
    [t5, y5]    = ode15s(@fn_Ecoli,tswof,y0,options,data);

    y0          = y5(end,:);                      
    Fe0         = fe;

    T           = [T;t4(2:end);t5(2:end)];         
    Y           = [Y;y4(2:end,:);y5(2:end,:)];
end
    %% Visualize profiles
      
Plotprofiles





