function int_msetupID = getMeasuringSetupID(int_runID, str_variableType, conn)
    qlquery = ['SELECT measuring_setup_id FROM measuring_setup WHERE run_id = ' num2str(int_runID) ' AND variable_type_id = (SELECT variable_type_id FROM variable_types WHERE canonical_name = ''' str_variableType ''');'];
    int_msetupID = SQLquery(qlquery, conn);
    int_msetupID = int_msetupID(1);
end



