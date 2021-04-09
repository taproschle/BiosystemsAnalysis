function int_currentCulivationAge = getCCultivationAge(int_runID, conn)
    str_sql = ['SELECT start_time FROM runs WHERE run_id = ' num2str(int_runID) ';'];
    
    int_currentCulivationAge = SQLquery(str_sql, conn);
end