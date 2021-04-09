%% Set preferences with setdbprefs.
setdbprefs('DataReturnFormat', 'cellarray');
setdbprefs('NullNumberRead', 'NaN');
setdbprefs('NullStringRead', 'null');
run_id = 75;

%% Make connection to database using JDBC driver and iLab credentials.

conn = database('ilabdb', 'Autobio', 'bvt1autobio!', 'Vendor', 'MYSQL', 'Server', 'ht-server.bioprocess.tu-berlin.de', 'PortNumber', 3306);

%% Test connection, read something from db
% sqlquery = 'Select * From ilabdb.setpoints where profile_id = 564';
% sqlquery = 'Select * From ilabdb.measurements_experiments Where experiment_id = 6757 and measuring_setup_id = 879';

sqlquery = 'Select * From ilabdb.runs';
curs = exec(conn,sqlquery); 
curs = fetch(curs);
All_runs = curs.Data

sqlquery1 = ['Select * From ilabdb.bioreactors Where run_id = ' num2str(run_id)];
curs = exec(conn,sqlquery1); 
curs = fetch(curs);
bioreactor_id = curs.Data(1,1);

sqlquery2 = ['Select * From ilabdb.experiments Where bioreactor_id = ' num2str(bioreactor_id)];
curs = exec(conn,sqlquery2); 
curs = fetch(curs);
profile_ids = curs.Data(:,4);


%% Use standard sql query language to write setpoints to ilab
sqlquery = 'Insert into ilabdb.setpoints (profile_id, variable_type_id, scope, cultivation_age, setpoint_value, checksum) Values (564, 25, ''e'', 3800, 220, 1)';
curs = exec(conn,sqlquery); % message displays error of child-parent-table foreign key constraint.
curs = fetch(curs);
curs.Data

%% Use 'datainsert' to write setpoints to ilab
colnames = {'profile_id', 'variable_type_id', 'scope', 'cultivation_age', 'setpoint_value', 'checksum'};

data = {564, 25, 'e', 3800, 220, 1};

% data_table = cell2table(data,'VariableNames',colnames)

tablename = 'ilabdb.setpoints';

datainsert(conn,tablename,colnames,data) % message displays error of child-parent foreign key constraint.


%{{{select child.id from child left join parent on (child.parent_id=parent.id) where child.id is not null and parent.id is null;}}}


%% Close database connection.
close(conn);
clear curs conn

