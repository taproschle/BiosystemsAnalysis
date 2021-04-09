function  Exp_Const_Pulse_Feed

%{
by Emmanuel Anane
TU-Berlin
Inst. of Biotechnology

2018-02-10

General notes:

 Pulse fed-batch feed calculation for robot cultivations
 No reconciliation of pH control acid/base additions with the total
 volume. acid/base volumes are assumed to compensate for evaporation.
 added time vector to output arguments--20170910
 added sampling effects in the constant feed phase--20170820
 mechanistic model: qO is limited by oxygen availability andd Ko
 mechanistic model: using the steady state DOT concept (reduced DOT from
 ODE system to simple algebraic system).

Part I in this code is only for coonecting to the iLab database on 
Emmanuel's mac environment, if you use other systems (Windows), you have
to replace this part with your own code that gives the OD600 vector and the 
cultivation_age_sec.


##########################################################################


 PART I:
 
- connect to ilab database
- read last measured OD600 values, put these in a vector called 'OD600' and 
  make sure length(OD600) = number_reactors
- read the cultivation age from ilab in hours.

%}
%% RE-EXECUTING CODE DURING A CULTIVATION
%{
Check if this is the first time the code is being run
during a cultivation. In subsequent executions of the code (e.g. after each
sampling point when new OD values are available), pulse setpoints that are
not yet executed by the robot are deleted from the database, and new ones
(obtained from calculations using updated analysis results) are written to 
the db.
%}
% prompt1         = 'Re-execution of the code?'; %input yes or no
% prompt2         = 'How many hours of exponential feed are left?';
% prompt3         = 'How many hours of constant feed are left?';
% 
% repeating_code          = input(prompt1,'s');
% Specs.exponential_time  = input(prompt2);              %out of the fed-batch time, how long is exponential feeding, in hours? The rest of the time is constant feeding.
% Specs.const_time        = input(prompt3);              %out of the fed-batch time, how long is exponential feeding, in hours? The rest of the time is constant feeding.

repeating_code          = 'no';
Specs.exponential_time  =  3;              %out of the fed-batch time, how long is exponential feeding, in hours? The rest of the time is constant feeding.
Specs.const_time        =  6;              %out of the fed-batch time, how long is exponential feeding, in hours? The rest of the time is constant feeding.


run_id          = 275;      % Get this from the 2mag interface on the robot


%% Connect to iLab database

% Set database preferences with setdbprefs.
setdbprefs('DataReturnFormat', 'cellarray'); %cell array imports everything in 'cell', and converting to 'double' is tedious, one can use cell2mat to convert.
setdbprefs('NullNumberRead', 'NaN');
setdbprefs('NullStringRead', 'null');


% Make connection to database using JDBC driver and iLab credentials.

conn = database('ilabdb', 'Autobio', 'bvt1autobio!', 'Vendor', 'MYSQL', 'Server', 'ht-server.bioprocess.tu-berlin.de', 'PortNumber', 3306);

%% Extract recently measured OD600 values from ilab database

sqlquery1 = ['SELECT bioreactor_id FROM ilabdb.bioreactors WHERE run_id = ' num2str(run_id)];
curs = exec(conn,sqlquery1); curs = fetch(curs);
bioreactor_id = cell2mat(curs.Data(1,1));

sqlquery2 = ['SELECT experiment_id, profile_id FROM ilabdb.experiments WHERE bioreactor_id = ' num2str(bioreactor_id) ' ORDER BY container_number;'];
curs = exec(conn,sqlquery2); curs = fetch(curs);
experiment_ids = cell2mat(curs.Data(:,1)); % for fetching OD600 values
profile_ids = cell2mat(curs.Data(:,2));% for writing setpoints

OD600  = zeros(1,length(experiment_ids));
int_MSID_OD600 = cell2mat(getMeasuringSetupID(run_id, 'OD600', conn)); % get the measurement id for OD600 values

for k0 = 1:length(experiment_ids)
    sqlquery3 = ['SELECT measured_value FROM ilabdb.measurements_experiments WHERE experiment_id = ' num2str(experiment_ids(k0)) ' AND measuring_setup_id = ' num2str(int_MSID_OD600)]; 
    curs = exec(conn,sqlquery3); curs = fetch(curs);
    OD600(k0) = cell2mat(curs.Data(end));
end
od1 = OD600(1,2:8); od2 = OD600(1,10:16);od3 = OD600(1,18:24);
ODnew = [od1,od2,od3];
Avg_OD = mean(ODnew); % use the mean OD to calculate the feed rates, thus forcing the all cultivations to start from the same point, and to ensure that the same amount of glucose is fed to each.
%% Extract cultivation time/age

start_time = getCCultivationAge(run_id, conn);
t1 = datevec(start_time);
t2 = clock;

cultivation_age = etime(t2,t1);    % etime calculates elapsed time between t1 and t2 in hours
cultivation_age = cultivation_age./3600;
% Example of OD600 and cultivation age as outputs from Part I,
OD600 = ones(1,length(profile_ids)).*Avg_OD;
% OD600 = [6.0 6.0 6.5 6.9 6.2 5.6 5.7 6.35];     %OD600 values at the end of batch phase.
OD600 = [6.92 6.92 6.88 6.88 7.36 7.04 6.68 6 6.8 6.8 7.08 7.2 7.12 7.08 7.08 7.4 7.2 7.2 7.32 6.88 6.84 6.88 6 6.12];

cultivation_age = 5.55;                        %cultivation age in hours

%{
###########################################################################

PART II

Actual calculation of pulse volumes. Give the specifications of your cultivation 
in the structure called 'Specs'.  

%}

%% User inputs for pulse feed
Specs.Si                    = 200;            %glucose concentration in feed solution, g/L
Specs.fed_batch_duration    = Specs.exponential_time +  Specs.const_time  ;          
Specs.mu_set                = 0.350;          % the set specific growth rate in the exponential feeding phase
Specs.T_vector              = [1 4];          % feeding (1min), resting time (9min)] total pulse cycle = 10min
Specs.pulses                = Specs.fed_batch_duration*60/sum(Specs.T_vector);  %total number of pulses to be given
Specs.V_reactor             = 10.*ones(1,length(OD600));    % reactor volumes in ml
Specs.V_sample              = 0.3*ones(1,length(OD600));   % volume of samples taken in ml 
Specs.Sample_times          = 0:12:Specs.pulses;             % A sample is taken after every 6th pulse == 1hr
Specs.plot_feed_profile     = 'yes';                        % put 'no' if you dont want a plot of the feed profiles.
Specs.double_pulses_int     = 'yes';                        % specify if some pulses are double spaced
Specs.Pulse_interval        = 2;
Specs.which_bioreactors_dbp = [6 7 8 14 15 16 22 23 24];
Specs.enzyme_addition       = 'yes';                       
Specs.which_bioreactors_enz = [2 10 18];
Specs.enz_Vfactor           = 5;                        % the factor by which glucose pulse should be divided to yield the enzyme volume, the enzyme is not being consumed, so additions should take into account what was added previously.
Specs.no_pulses             = 'yes';                       
Specs.which_bioreactors_np  = [1 9 17];
All_bioreactors             = 1:1:length(OD600);


[Pulse_time,V_pipette, Total_vol_added]    = Calculate_pulses(Specs,OD600,cultivation_age);

disp(Total_vol_added)

% enz_pulses = Specs.exponential_time*60/sum(Specs.T_vector)

%{
###########################################################################

PART III: Write feed volumes as setpoints to ilab

The 'datainsert' function of matlab is used here. Alternatively, sql 'insert
into' query can also be used to write the setpoints to ilab. If some
bioreactors are excluded from the pulse feeding, the indices k3 and k4 can
be changed accordingly.

%}

if strcmp(repeating_code,'yes') % delete future setpoints of previous code-execution
    for k2 = 1:length(OD600)
       sqlquery =  ['DELETE FROM ilabdb.setpoints WHERE '... 
                    'profile_id = ' num2str(profile_ids(k2))...
                    ' AND cultivation_age >= ' num2str(Pulse_time(1))...
                    (' AND variable_type_id = 87 OR variable_type_id = 70')];
       exec(conn,sqlquery);
    end
end   


colnames = {'profile_id', 'variable_type_id', 'scope', 'cultivation_age', 'setpoint_value', 'checksum'};

%----------------- write glucose setpoints--------------------------------
glu_pulse_rxtors = setdiff(All_bioreactors,Specs.which_bioreactors_enz);

profile_id_glu = profile_ids(glu_pulse_rxtors);

    for k3 = 1:length(glu_pulse_rxtors)

        profile_id = profile_id_glu(k3);

        for k4 = 1:length(V_pipette)

            feed_data = {profile_id, 87, 'e', Pulse_time(k4), V_pipette(k4,glu_pulse_rxtors(k3)), 1};

            tablename = 'ilabdb.setpoints';

            datainsert(conn,tablename,colnames,feed_data)
        end
    end
    
%----------------- write enzyme setpoints--------------------------------
if strcmp(Specs.enzyme_addition,'yes') % && strcmp(repeating_code,'no')
    
    profile_ids_enz = profile_ids(Specs.which_bioreactors_enz);
    
    for k5 = 1:length(profile_ids_enz)

        profile_id = profile_ids_enz(k5);

        for k6 = 1 : Specs.exponential_time*60/sum(Specs.T_vector) % give enzyme pulses for only exponential feeding phase

            feed_data = {profile_id, 70, 'e', Pulse_time(k6), V_pipette(k6,Specs.which_bioreactors_enz(k5)), 1};

            tablename = 'ilabdb.setpoints';

            datainsert(conn,tablename,colnames,feed_data)
        end
    end
end
 close(conn)
end

% Hurray!!!!  Model-based feed profiles.




