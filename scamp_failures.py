#!/usr/bin/env python

import coreutils.desdbi


def connectDB(dbase):
    """
    Connect to database, query and extract data from table, and close connection to database
    dbase can be: 
    db-destest
    db-desoper
    """
    try:
        desdmfile = os.environ["des_service"]
    except KeyError:
        desdmfile  = None

    dbh = coreutils.desdbi.DesDbi(desdmfile,dbase)
    #cur = dbh.cursor()
        
    return dbh

def query_to_cur(dbh, query, verbose=verbose):
    """Execute a query and return a cursor to a query
    :param dbh: Database connection handler
    :param query: query to execute
    :param debug: verbosity
    
    """
    if verbose : print query
    cur = dbh.cursor()
    cur.execute(query)
    return cur


def main():
    

    section = 'db-desoper'
    dbh = connectDB(section)

    query = """ SELECT distinct exposure.exposurename,exposure.object,image.band,image.exposureid,location.nite,location.run,
                run_data_state.state FROM exposure,location,image,run_data_state WHERE location.fileclass='red' AND 
                location.filetype='red' AND location.project='OPS' AND REGEXP_LIKE(location.archivesites,'[^N]') AND 
                location.id=image.id AND image.scampflg='1' AND run_data_state.run=location.run and run_data_state.state='NEW' 
                and image.exposureid=exposure.id  and exposure.object like '\%SN-\%'  order by location.nite;"""
    results = query_to_cur(dbh,query)
    
    for values in results:
        print values

    return 0

if  __name__ == '__main__':
    status = main()
    sys.exit(status)


