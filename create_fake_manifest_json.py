#!/usr/bin/env python

import coreutils.desdbi
import os

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
    if dbh.is_postgres():
        print 'Connected to postgres DB'
    elif dbh.is_oracle():
        print 'Connected to oracle DB'
        print 'which_services_file = ', dbh.which_services_file()
        print 'which_services_section = ', dbh.which_services_section()
        
    cur = dbh.cursor()
        
    return cur

def main():
    
    
    query_items = ['expnum','nite', 'band', 'object', 'propid', 'date_obs']
    
    data = {}
    for index, item in enumerate(query_items):
        data[item] = index
    
    print data
    
    query = "select expnum, nite, band, object, propid, date_obs from exposure where nite='20121204' and object like '%%SN-%%'"
    cur = connectDB("db-desoper")
    
    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)

    expnum = []
    nite = []
    band = []
    object = []
    propid = []
    date_obs = []
    for item in cur:
        print item
        expnum=data['expnum'] = item[0]
        data['nite'] = item[1]
        data['band'] = item[2]
        data['object'] = item[3]
        data['propid'] = item[4]
        data['date_obs']= item[5]

    print data
    
    
    
if  __name__ == '__main__':
    main()
           