function str_query = SQLquery(str_sql, conn)
    curs = exec(conn,str_sql); 
    curs = fetch(curs);
    str_query = curs.Data;
end