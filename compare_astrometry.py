#!/usr/bin/env python
"""
Plot different parameters stored on ricardoc.flats_stats_qa

.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>

  Needs to setup the following before running:
  matplotlib
  
"""
# imports
#import coreutils.desdbi
import argparse
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import datetime as dt
import matplotlib.dates as mdates
import numpy as np
import os
import coreutils.desdbi
import cx_Oracle
import sys


parser = argparse.ArgumentParser(description="Compare astrometry from 2013 and reprocessing 2014 for SNe fileds")
parser.add_argument("-v", "--verbose", required=False,action="store_true",default=False,help="Use --verbose for verbose output")

args = parser.parse_args()
verbose = args.verbose


def connectDB(dbase):
    """
    Connect to database, query and extract data from table, and close connection to database
    dbase can be: 
    db-destest
    db-desoper
    """
    try:
        desdmfile = os.environ["DES_SERVICES"]
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
    if verbose : 
        print query
    cur = dbh.cursor()
    #cur.array
    cur.execute(query)
    #print "Executing query :", query_flats
    #dbcursor.arraysize = 1000 # get 1000 at a time when fetching
    #dbcursor.execute(query_flats)

    return cur


def main():



    section = 'db-desoper'
    dbh = connectDB(section)



    query_Y1A1_20131112 = """select distinct(exposureid),run,nite,ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                from exposure_qa where run = '20140513162817_20131112' """

    query_snese_20131112 = """select distinct(exposureid),run,nite,ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                from exposure_qa where run = '20140510135337_20131112' """
    
    
    data_Y1A1 = query_to_cur(dbh, query_Y1A1_20131112, verbose=True)
    data_snese = query_to_cur(dbh, query_snese_20131112, verbose=True)


    expid_Y1A1_2014 = []
    run_Y1A1_2014 = []
    nite_Y1A1_2014 = []
    rms_1_Y1A1_2014 = []
    rms_2_Y1A1_2014 = []
    rms_hi_1_Y1A1_2014 = []
    rms_hi_2_Y1A1_2014 = []
    
    for val in data_Y1A1:
        expid_Y1A1_2014.append(val[0])
        run_Y1A1_2014.append(val[1])
        nite_Y1A1_2014.append(val[2])
        rms_1_Y1A1_2014.append(val[3])
        rms_2_Y1A1_2014.append(val[4])
        rms_hi_1_Y1A1_2014.append(val[5])
        rms_hi_2_Y1A1_2014.append(val[6])
        

    expid_snese_2014 = []
    run_snese_2014 = []
    nite_snese_2014 = []
    rms_1_snese_2014 = []
    rms_2_snese_2014 = []
    rms_hi_1_snese_2014 = []
    rms_hi_2_snese_2014 = []
    
    for val in data_snese:
        expid_snese_2014.append(val[0])
        run_snese_2014.append(val[1])
        nite_snese_2014.append(val[2])
        rms_1_snese_2014.append(val[3])
        rms_2_snese_2014.append(val[4])
        rms_hi_1_snese_2014.append(val[5])
        rms_hi_2_snese_2014.append(val[6])
        
        
    plt.clf()
    #plt.xticks(rotation=45)
    #plt.grid()
    plt.rc("font", size=8)

    plt.figure(1)
    plt.subplot(2,1,1)
    #plt.plot(rms_1_snese_2014,rms_1_Y1A1_2014,'bo')
    n, bins, patches = plt.hist(rms_1_snese_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', linestyle=('solid'),lw=2, color='blue', label='SNe 2014 params rms_1 low')
    n, bins, patches = plt.hist(rms_1_Y1A1_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', linestyle=('solid'),lw=2, color='red', label='Y1A1 params rms_1 low')
    plt.xlabel('ASTROMSIGMA_REF_1')
    plt.ylabel('Nun Exposures')
    plt.legend(loc='upper right')

    plt.subplot(2,1,2)
    n, bins, patches = plt.hist(rms_2_snese_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', linestyle=('solid'), lw=2, color='blue', label='SNe 2014 params rms_2 low')
    n, bins, patches = plt.hist(rms_2_Y1A1_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', linestyle=('solid'),lw=2, color='red', label='Y1A1 params rms_2 low')
    plt.xlabel('ASTROMSIGMA_REF_2')
    plt.ylabel('Nun Exposures')

    #n, bins, patches = plt.hist(rms_hi_1_snese_2014, bins=40, range=[0,1] ,histtype = 'step', lw=1, color='green', label='SNe 2014 params rms Hi')
    #n, bins, patches = plt.hist(rms_hi_1_Y1A1_2014, bins=40, range=[0,1] ,histtype = 'step', lw=1, color='black', label='Y1A1 params rms Hi')

    #plt.plot(rms_2_snese_2014,rms_2_Y1A1_2014,'ro')
    plt.legend(loc='upper right')
    plt.savefig('Y1A1_params_vs_SNe_params_nov12.pdf')
    plt.close()

    plt.show()


    plt.figure(1)
    plt.subplot(2,1,1)
    n, bins, patches = plt.hist(rms_1_snese_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', lw=2, color='blue', label='SNe 2014 params rms_1 Hi')
    n, bins, patches = plt.hist(rms_1_Y1A1_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', lw=2, color='red', label='Y1A1 params rms_1 Hi')
    plt.xlabel('ASTROMSIGMA_REF_HIGHSN_1')
    plt.ylabel('Num Exposures')
    plt.legend(loc='upper right')

    
    plt.subplot(2,1,2)
    n, bins, patches = plt.hist(rms_2_snese_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', lw=2, color='blue', label='SNe 2014 params rms_2 Hi')
    n, bins, patches = plt.hist(rms_2_Y1A1_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', lw=2, color='red', label='Y1A1 params rms_2 Hi')
    plt.xlabel('ASTROMSIGMA_REF_HIGHSN_1')
    plt.ylabel('Num Exposures')
    plt.legend(loc='upper right')
    
    
    plt.show()

    """
    query_Margaret_test = select ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2 from exposure_qa where
                            run = '20140513151817_20130902' and  exposureid in (select distinct(exposureid) from image where run = '20140513151817_20130902' and project = 'ACT' 
                            and exposureid in (select id from exposure where object like '%%SN-%%' and nite='20130902'))
    
    
    data_Margaret = query_to_cur(dbh, query_Margaret_test, verbose=True)

    rms_1_test_2014 = []
    rms_2_test_2014 = []
    rms_hi_1_test_2014 = []
    rms_hi_2_test_2014 = []
    
    for val in data_Margaret:
        rms_1_test_2014.append(val[0])
        rms_2_test_2014.append(val[1])
        rms_hi_1_test_2014.append(val[2])
        rms_hi_2_test_2014.append(val[3])

    plt.figure(1)
    plt.subplot(2,1,1)
    n, bins, patches = plt.hist(rms_1_snese_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', lw=2, color='blue', label='SNe 2014  rms_1')
    n, bins, patches = plt.hist(rms_1_test_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', lw=2, color='red', label='Y1A1 TEST params rms_1')
    plt.xlabel('ASTROMSIGMA_REF_1')
    plt.ylabel('Num Exposures')
    plt.legend(loc='upper right')

    plt.subplot(2,1,2)
    n, bins, patches = plt.hist(rms_2_snese_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', lw=2, color='blue', label='SNe 2014 params rms_2')
    n, bins, patches = plt.hist(rms_2_test_2014, bins=40, range=[0.1,0.6] ,histtype = 'step', lw=2, color='red', label='Y1A1 TEST params rms_2')
    plt.xlabel('ASTROMSIGMA_REF_2')
    plt.ylabel('Num Exposures')
    plt.legend(loc='upper right')

    """
    
    #Another night to check in the morning:
    query_sept7_margaret = """select distinct(exposureid),ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2 from exposure_qa where run='20140513194115_20130907'  
                                and  exposureid in (select id from exposure where object like '%%SN-%%' and id in (select exposureid from image where
                                run='20140513194115_20130907' and project='ACT')) order by exposureid"""

    query_sept2_margaret = """select distinct(exposureid),ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2 from exposure_qa where run = '20140513151817_20130902'  
                                and  exposureid in (select id from exposure where object like '%%SN-%%' and id in (select exposureid from image where
                                run = '20140513151817_20130902' and project='ACT')) order by exposureid"""
    

    query_sept11_margaret = """select distinct(exposureid),ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2 from exposure_qa where run = '20140514090259_20130911'  
                                and  exposureid in (select id from exposure where object like '%%SN-%%' and id in (select exposureid from image where
                                run = '20140514090259_20130911' and project='ACT')) order by exposureid"""

    
    #same thing with sept 09 20140513194439_20130909
    query_sep7_snese_2014 = """select distinct(exposureid),ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                from exposure_qa where run = '20140425140807_20130907'  order by exposureid"""

    query_sep2_snese_2014 = """select distinct(exposureid),ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                from exposure_qa where run='20140425135930_20130902'  order by exposureid"""

    query_sep11_snese_2014 = """select distinct(exposureid),ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                from exposure_qa where run='20140429112229_20130911'  order by exposureid"""

    
    data_y1a1_sep7 = query_to_cur(dbh,query_sept7_margaret , verbose=True)
    data_snese_sep2 = query_to_cur(dbh, query_sep2_snese_2014, verbose=True)
    data_y1a1_sep2 = query_to_cur(dbh,query_sept2_margaret , verbose=True)
    data_snese_sep7 = query_to_cur(dbh, query_sep7_snese_2014, verbose=True)
    data_y1a1_sep11 = query_to_cur(dbh,query_sept11_margaret , verbose=True)
    data_snese_sep11 = query_to_cur(dbh, query_sep11_snese_2014, verbose=True)


    expid_test_sep7_2014 = []
    rms_1_test_sep7_2014 = []
    rms_2_test_sep7_2014 = []
    rms_hi_1_test_sep7_2014 = []
    rms_hi_2_test_sep7_2014 = []
    
    for val in data_y1a1_sep7:
        expid_test_sep7_2014.append(val[0])
        rms_1_test_sep7_2014.append(val[1])
        rms_2_test_sep7_2014.append(val[2])
        rms_hi_1_test_sep7_2014.append(val[3])
        rms_hi_2_test_sep7_2014.append(val[4])

    expid_test_sep2_2014 = []
    rms_1_test_sep2_2014 = []
    rms_2_test_sep2_2014 = []
    rms_hi_1_test_sep2_2014 = []
    rms_hi_2_test_sep2_2014 = []
    

    for val in data_y1a1_sep2:
        expid_test_sep2_2014.append(val[0])
        rms_1_test_sep2_2014.append(val[1])
        rms_2_test_sep2_2014.append(val[2])
        rms_hi_1_test_sep2_2014.append(val[3])
        rms_hi_2_test_sep2_2014.append(val[4])


    expid_test_sep11_2014 = []
    rms_1_test_sep11_2014 = []
    rms_2_test_sep11_2014 = []
    rms_hi_1_test_sep11_2014 = []
    rms_hi_2_test_sep11_2014 = []
    

    for val in data_y1a1_sep11:
        expid_test_sep11_2014.append(val[0])
        rms_1_test_sep11_2014.append(val[1])
        rms_2_test_sep11_2014.append(val[2])
        rms_hi_1_test_sep11_2014.append(val[3])
        rms_hi_2_test_sep11_2014.append(val[4])


    expid_snese_2014 = []
    rms_1_snese_2014 = []
    rms_2_snese_2014 = []
    rms_hi_1_snese_2014 = []
    rms_hi_2_snese_2014 = []
    
    for val in data_snese_sep7:
        expid_snese_2014.append(val[0])
        rms_1_snese_2014.append(val[1])
        rms_2_snese_2014.append(val[2])
        rms_hi_1_snese_2014.append(val[3])
        rms_hi_2_snese_2014.append(val[4])
   
    #expid_snese_2014 = []
    #rms_1_snese_2014 = []
    #rms_2_snese_2014 = []
    #rms_hi_1_snese_2014 = []
    #rms_hi_2_snese_2014 = []
   
   
    for val in data_snese_sep2:
        expid_snese_2014.append(val[0])
        rms_1_snese_2014.append(val[1])
        rms_2_snese_2014.append(val[2])
        rms_hi_1_snese_2014.append(val[3])
        rms_hi_2_snese_2014.append(val[4])

    for val in data_snese_sep11:
        expid_snese_2014.append(val[0])
        rms_1_snese_2014.append(val[1])
        rms_2_snese_2014.append(val[2])
        rms_hi_1_snese_2014.append(val[3])
        rms_hi_2_snese_2014.append(val[4])

    
    
    #Select only the rms from the exposures from 0902 and 0907 reduced in snese
    allids_y1a1_sep7 = ','.join(str(id) for id in expid_test_sep7_2014)
    allids_y1a1_sep2 = ','.join(str(id) for id in expid_test_sep2_2014)
    allids_y1a1_sep11 = ','.join(str(id) for id in expid_test_sep11_2014)

    query_snese_7_2014_ricardo = """select exposureid,ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                                from exposure_qa where run = '20140425140807_20130907' and exposureid in (%s)""" % allids_y1a1_sep7
    
    query_snese_2_2014_ricardo = """select exposureid,ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                                from exposure_qa where run = '20140425135930_20130902' and exposureid in (%s)""" % allids_y1a1_sep2

    query_snese_11_2014_ricardo = """select exposureid,ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                                from exposure_qa where run = '20140429112229_20130911' and exposureid in (%s)""" % allids_y1a1_sep11


    data_snese_matched_sep7 = query_to_cur(dbh,query_snese_7_2014_ricardo , verbose=True)
    data_snese_matched_sep2 = query_to_cur(dbh,query_snese_2_2014_ricardo, verbose=True)
    data_snese_matched_sep11 = query_to_cur(dbh,query_snese_11_2014_ricardo, verbose=True)

    expid_snese_matched_2014 = []
    rms_1_snese_matched_2014 = []
    rms_2_snese_matched_2014 = []
    rms_hi_1_snese_matched_2014 = []
    rms_hi_2_snese_matched_2014 = []
    
    for val in data_snese_matched_sep7:
        expid_snese_matched_2014.append(val[0])
        rms_1_snese_matched_2014.append(val[1])
        rms_2_snese_matched_2014.append(val[2])
        rms_hi_1_snese_matched_2014.append(val[3])
        rms_hi_2_snese_matched_2014.append(val[4])

    for val in data_snese_matched_sep2:
        expid_snese_matched_2014.append(val[0])
        rms_1_snese_matched_2014.append(val[1])
        rms_2_snese_matched_2014.append(val[2])
        rms_hi_1_snese_matched_2014.append(val[3])
        rms_hi_2_snese_matched_2014.append(val[4])

    for val in data_snese_matched_sep11:
        expid_snese_matched_2014.append(val[0])
        rms_1_snese_matched_2014.append(val[1])
        rms_2_snese_matched_2014.append(val[2])
        rms_hi_1_snese_matched_2014.append(val[3])
        rms_hi_2_snese_matched_2014.append(val[4])

    
    #print rms_1_test_2014,rms_1_test_2014,rms_1_snese_2014,rms_2_snese_2014

    plt.figure(1)
    plt.subplot(2,1,1)
    n, bins, patches = plt.hist(rms_1_snese_2014, bins=40, range=[0,0.6] ,histtype = 'step', lw=2, color='blue', label='SNeSE 2014 0902,0907,0911 rms_1')
    n, bins, patches = plt.hist(rms_1_test_sep2_2014+rms_1_test_sep7_2014+rms_1_test_sep11_2014, bins=40, range=[0,0.6] ,histtype = 'step', lw=2, color='red', label='Y1A1 Test 0902,0907,0911 rms_1')
    #n, bins, patches = plt.hist(rms_1_test_sep7_2014, bins=40, range=[0,0.6] ,histtype = 'step', lw=2, color='red', label='Y1A1 Test 0902,0907,0911 rms_1')
    #n, bins, patches = plt.hist(rms_1_test_sep11_2014, bins=40, range=[0,0.6] ,histtype = 'step', lw=2, color='red', label='Y1A1 Test 0902,0907,0911 rms_1')
    n, bins, patches = plt.hist(rms_1_snese_matched_2014, bins=40, range=[0,0.6] ,histtype = 'step', lw=2, color='green',linestyle='dashed', label='SNeSE 0902,0907,0911 rms_1')
    plt.xlabel('ASTROMSIGMA_REF_1')
    plt.ylabel('Num Exposures')
    plt.title('nights 0902-0911')
    plt.ylim(0,20)
    plt.legend(loc='upper right')

    
    plt.subplot(2,1,2)
    n, bins, patches = plt.hist(rms_2_snese_2014, bins=40, range=[0,0.6] ,histtype = 'step', lw=2, color='blue', label='SNeSE 2014 0902,0907,0911 rms_2')
    n, bins, patches = plt.hist(rms_2_test_sep2_2014+rms_2_test_sep7_2014+rms_2_test_sep11_2014, bins=40, range=[0,0.6] ,histtype = 'step', lw=2, color='red', label='Y1A1 Test 0902,0907,0911 rms_2')
    #n, bins, patches = plt.hist(rms_2_test_sep7_2014, bins=40, range=[0,0.6] ,histtype = 'step', lw=2, color='red', label='Y1A1 Test 0902,0907,0911 rms_2')
    #n, bins, patches = plt.hist(rms_2_test_sep11_2014, bins=40, range=[0,0.6] ,histtype = 'step', lw=2, color='red', label='Y1A1 Test 0902,0907,0911 rms_2')
    n, bins, patches = plt.hist(rms_2_snese_matched_2014, bins=40, range=[0,0.6] ,histtype = 'step', lw=2, color='green',linestyle='dashed', label='SNeSE 0902,0907,0911 rms_2')
    plt.xlabel('ASTROMSIGMA_REF_2')
    plt.ylabel('Num Exposures')
    plt.title('nights 0902-0911')
    plt.ylim(0,20)
    plt.legend(loc='upper right')

    #plt.show()
    plt.savefig('Y1A1_latest_vs_SNe_processing.pdf')
    plt.close()

    
    plt.show()

    allexps_Y1A1 = expid_test_sep2_2014 + expid_test_sep7_2014 + expid_test_sep11_2014
    all_rms_1_test = rms_1_test_sep2_2014 + rms_1_test_sep7_2014 + rms_1_test_sep11_2014

    all_rms1_exps = []
    for i, val in enumerate(all_rms_1_test):
        if val > 0.25:
            all_rms1_exps.append(allexps_Y1A1[i])
            
    expids = ','.join(str(id) for id in all_rms1_exps)
    
    query_allexps = """select id,nite,exposurename,band,exptime,object from exposure where id in (%s)""" % expids
    data_allexps = query_to_cur(dbh,query_allexps, verbose=True)

    for val in data_allexps:
        print val

    
    #allexpids_margaret = ','.join(str(id) for id in expid_test_2014)
    #allexpids_ricardo = ','.join(str(id) for id in expid_snese_2014)
    #query_margaret = """ select exposurename, nite,band,object from exposure where id in (%s)""" % allexpids_margaret
    #query_ricardo = """ select exposurename, nite,band,object from exposure where id in (%s)""" % allexpids_ricardo

    #data_Y1A1 = query_to_cur(dbh,query_margaret , verbose=True)
    #data_ricardo = query_to_cur(dbh, query_ricardo, verbose=True)

    #for val in data_Y1A1:
    #    print 'Y1A1...', val
    #for val in data_ricardo:
    #    print 'snese...', val
    
    

    query_2014 = """select distinct(exposureid),run,nite,ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                from exposure_qa where run in
                (select r.run from run r, run_data_state rt
                where r.pipeline='snese'
                and r.project='OPS'
                and r.run_submit > '24-APR-14'
                and r.run_submit < '11-MAY-14'
                and rt.run = r.run
                and rt.state != 'JUNK')
                order by nite"""

    
    data_2014 = query_to_cur(dbh, query_2014, verbose=False)


    expid_2014 = []
    run_2014 = []
    nite_2014 = []
    rms_1_2014 = []
    rms_2_2014 = []
    rms_hi_1_2014 = []
    rms_hi_2_2014 = []
    
    for val in data_2014:
        expid_2014.append(val[0])
        run_2014.append(val[1])
        nite_2014.append(val[2])
        rms_1_2014.append(val[3])
        rms_2_2014.append(val[4])
        rms_hi_1_2014.append(val[5])
        rms_hi_2_2014.append(val[6])
        
    print len(expid_2014)
    if len(expid_2014) >= 900:
        set_1 = expid_2014[0:900]
        set_2 = expid_2014[900:]
        
    allexps_set1 = ','.join(str(x) for x in set_1)
    allexps_set2 = ','.join(str(x) for x in set_2)

    
    query_2013_set1 = """select distinct(exposureid),run,nite,ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                from exposure_qa where run in
                (select r.run from run r, run_data_state rt
                where r.pipeline='snese'
                and r.project='OPS'
                and r.run_submit > '31-AUG-13'
                and r.run_submit < '13-NOV-13'
                and rt.run = r.run
                and rt.state != 'JUNK')
                and exposureid in (%s)
                order by nite""" % (allexps_set1)

    query_2013_set2 = """select distinct(exposureid),run,nite,ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                from exposure_qa where run in
                (select r.run from run r, run_data_state rt
                where r.pipeline='snese'
                and r.project='OPS'
                and r.run_submit > '31-AUG-13'
                and r.run_submit < '13-NOV-13'
                and rt.run = r.run
                and rt.state != 'JUNK')
                and exposureid in (%s)
                order by nite""" % (allexps_set2)


    data_2013_set1 = query_to_cur(dbh, query_2013_set1, verbose=False)
    data_2013_set2 = query_to_cur(dbh, query_2013_set1, verbose=False)



    expid_2013 = []
    run_2013 = []
    nite_2013 = []
    rms_1_2013 = []
    rms_2_2013 = []
    rms_hi_1_2013 = []
    rms_hi_2_2013 = []
    
    for val in data_2013_set1:
        expid_2013.append(val[0])
        run_2013.append(val[1])
        nite_2013.append(val[2])
        rms_1_2013.append(val[3])
        rms_2_2013.append(val[4])
        rms_hi_1_2013.append(val[5])
        rms_hi_2_2013.append(val[6])

    for val in data_2013_set2:
        expid_2013.append(val[0])
        run_2013.append(val[1])
        nite_2013.append(val[2])
        rms_1_2013.append(val[3])
        rms_2_2013.append(val[4])
        rms_hi_1_2013.append(val[5])
        rms_hi_2_2013.append(val[6])



    query_Y1N_set1 = """select distinct(exposureid),run,nite,ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                from exposure_qa where run in
                (select r.run from run r, runtag rt
                where r.pipeline='firstcut'
                and r.project='OPS'
                and r.run_submit > '31-AUG-13'
                and r.run_submit < '13-NOV-13'
                and rt.run = r.run
                and rt.tag = 'Y1N_FIRSTCUT')
                and exposureid in (%s)
                order by nite""" % (allexps_set1)

    query_Y1N_set2 = """select distinct(exposureid),run,nite,ASTROMSIGMA_REF_1,ASTROMSIGMA_REF_2,ASTROMSIGMA_REF_HIGHSN_1,ASTROMSIGMA_REF_HIGHSN_2
                from exposure_qa where run in
                (select r.run from run r, runtag rt
                where r.pipeline='firstcut'
                and r.project='OPS'
                and r.run_submit > '31-AUG-13'
                and r.run_submit < '13-NOV-13'
                and rt.run = r.run
                and rt.tag = 'Y1N_FIRSTCUT')
                and exposureid in (%s)
                order by nite""" % (allexps_set2)

    
    
    data_Y1N_set1 = query_to_cur(dbh, query_Y1N_set1, verbose=False)
    data_Y1N_set2 = query_to_cur(dbh, query_Y1N_set2, verbose=False)

    expid_Y1N = []
    run_Y1N = []
    nite_Y1N = []
    rms_1_Y1N = []
    rms_2_Y1N = []
    rms_hi_1_Y1N = []
    rms_hi_2_Y1N = []
    
    for val in data_Y1N_set1:
        expid_Y1N.append(val[0])
        run_Y1N.append(val[1])
        nite_Y1N.append(val[2])
        rms_1_Y1N.append(val[3])
        rms_2_Y1N.append(val[4])
        rms_hi_1_Y1N.append(val[5])
        rms_hi_2_Y1N.append(val[6])

    for val in data_Y1N_set2:
        expid_Y1N.append(val[0])
        run_Y1N.append(val[1])
        nite_Y1N.append(val[2])
        rms_1_Y1N.append(val[3])
        rms_2_Y1N.append(val[4])
        rms_hi_1_Y1N.append(val[5])
        rms_hi_2_Y1N.append(val[6])



        
    #Plot a histogram with the astrometric soluctions form last year and this year.
    #Make sure I select the same exposures id so I can properly compare
    
    print "2013", len(expid_2014)
    print "2014", len(expid_2013)

    plt.figure(1)
    plt.subplot(2,1,1)

    #make a plot for each ccd        
    #plt.clf()
    #plt.xticks(rotation=45)
    #plt.grid()
    plt.rc("font", size=9)

    #print len(rms_1_Y1N)
    #print len(rms_1_2014)
    #n, bins, patches = plt.hist(rms_1_Y1N, bins=40, range=[0,1] ,histtype = 'step', lw=1, color='blue', label='Y1N rms low')
    n, bins, patches = plt.hist(rms_1_2014, bins=40, range=[0,0.6] ,histtype='step', lw=2, color='red', label='SNeSE Repro rms_1 low')
    n, bins, patches = plt.hist(rms_1_2013, bins=40, range=[0,0.6], histtype='step', lw=2, color='green', label='SNeSE 2013 rms_1 low')
    plt.xlabel('ASTROMSIGMA_REF_1')
    plt.ylabel('Num Exposures')
    #plt.title('Y1 vs Reprocessing')

    plt.legend(loc='upper right')


    plt.subplot(2,1,2)

    n, bins, patches = plt.hist(rms_2_2014, bins=40, range=[0,0.6] ,histtype='step', lw=2, color='red', label='SNeSE Repro rms_2 low')
    n, bins, patches = plt.hist(rms_2_2013, bins=40, range=[0,0.6], histtype='step', lw=2, color='green', label='SNeSE 2013 rms_2 low')
    plt.xlabel('ASTROMSIGMA_REF_2')
    plt.ylabel('Num Exposures')
    #plt.title('Y1 vs Reprocessing')
    plt.legend(loc='upper right')
    
    plt.show()
    
    plt.figure(1)
    plt.subplot(2,1,1)
    #n, bins, patches = plt.hist(rms_hi_1_Y1N,bins=40, range=[0,1], histtype = 'step', lw=1, color='blue', label='Y1N rms Hi')
    n, bins, patches = plt.hist(rms_hi_1_2014, bins=40, range=[0,0.6], histtype='step', lw=2, color='red', label='SNeSE Repro rms_1 Hi')
    n, bins, patches = plt.hist(rms_hi_1_2013, bins=40, range=[0,0.6], histtype='step', lw=2, color='green', label='SNeSE 2013 rms_1 Hi')
    plt.xlabel('ASTROMSIGMA_REF_HIGHSN_1')
    plt.ylabel('Num Exposures')
    #plt.title('Y1 vs Reprocessing')

    plt.legend(loc='upper right')

    
    plt.subplot(2,1,2)
    #n, bins, patches = plt.hist(rms_hi_1_Y1N,bins=40, range=[0,1], histtype = 'step', lw=1, color='blue', label='Y1N rms Hi')
    n, bins, patches = plt.hist(rms_hi_2_2014, bins=40, range=[0,0.6], histtype='step', lw=2, color='red', label='SNeSE Repro rms_2 Hi')
    n, bins, patches = plt.hist(rms_hi_2_2013, bins=40, range=[0,0.6], histtype='step', lw=2, color='green', label='SNeSE 2013 rms_2 Hi')
    plt.xlabel('ASTROMSIGMA_REF_HIGHSN_2')
    plt.ylabel('Num Exposures')
    #plt.title('Y1 vs Reprocessing')

    plt.legend(loc='upper right')

    
    
    plt.show()
    
    
    
    sys.exit()
    
    rms_lims = [-5,15]
    med_lims = [-5,15]
    if ccdnum in bkp5:
        rms_lims = [-5,15]
        med_lims = [-7,15]
        
    if ccdnum == 12 or ccdnum == 45:
        rms_lims = [0,30]
        
    #continue    
    fig=plt.figure()
    for i in range(1,7):
        ax = fig.add_subplot(3,2,i)
        if i == 1:
            #for a,ar,aw,ra,rar,raw in zip(med_a,med_ar,med_aw, rms_a, rms_ar,rms_aw):
            ax.plot(np.asarray(med_a), np.asarray(med_a) + np.asarray(rms_a), 'bo')
            ax.set_xlabel('med a',fontsize=9)
            ax.set_ylabel('med a + rms_a',fontsize=9)
            ax.set_title('CCD = %i' % ccdnum, fontsize=9)
            if ccdnum == 2 or ccdnum == 7:
                thetitle = 'CCD = %i, RMS_A= %3.1f, RMS_B=%3.1f' % (ccdnum, medianRMS_A, medianRMS_B)
                ax.set_title(thetitle, fontsize=7)
            ax.set_ylim(rms_lims)
            ax.set_xlim(med_lims)
            plt.axhline(y = cut_a, color = '#299967', ls = '--')
            for i in indx_a:
                ax.plot(np.asarray(med_a[i]), np.asarray(rms_a[i])+np.asarray(med_a[i]), 'r^', markersize=6)
            for i in indx_b:
                ax.plot(np.asarray(med_a[i]), np.asarray(rms_a[i])+np.asarray(med_a[i]), 'w^', markersize=6)
            #print "done plot 1"
        elif i == 2:
            ax.plot(np.asarray(med_ar), np.asarray(rms_ar) + np.asarray(med_ar), 'bo')
            ax.set_xlabel('med ar',fontsize=9)
            ax.set_ylabel('med ar + rms ar',fontsize=9)
            ax.set_ylim(rms_lims)
            plt.axhline(y = cut_a, color = '#299967', ls = '--')
            ax.set_xlim(med_lims)
            #print "done plot 2"
        elif i == 3:
            ax.plot(np.asarray(med_aw), np.asarray(rms_aw) + np.asarray(med_aw), 'bo')
            ax.set_xlabel('med aw',fontsize=9)
            ax.set_ylabel('med aw + rms aw',fontsize=9)
            ax.set_ylim(rms_lims)
            plt.axhline(y = cut_a, color = '#299967', ls = '--')
            ax.set_xlim(med_lims)
            #print "done plot 3"    
        elif i == 4:
            #for a,ar,aw,ra,rar,raw in zip(med_a,med_ar,med_aw, rms_a, rms_ar,rms_aw):
            ax.plot(np.asarray(med_b), np.asarray(rms_b) + np.asarray(med_b), 'bo')
            ax.set_xlabel('med b',fontsize=9)
            ax.set_ylabel('med b + rms b',fontsize=9)
            ax.set_ylim(rms_lims)
            plt.axhline(y = cut_b, color = '#299967', ls = '--')
            ax.set_xlim(med_lims)
            for i in indx_a:
                ax.plot(np.asarray(med_b[i]), np.asarray(rms_b[i])+np.asarray(med_b[i]), 'r^')
            for i in indx_b:
                ax.plot(np.asarray(med_b[i]), np.asarray(rms_b[i])+np.asarray(med_b[i]), 'w^',  markersize=5)
            #print "done plot 4"
        elif i == 5:
            ax.plot(np.asarray(med_br), np.asarray(rms_br) + np.asarray(med_br), 'bo')
            ax.set_xlabel('med br',fontsize=9)
            ax.set_ylabel('rms',fontsize=9)
            ax.set_ylim(rms_lims)
            plt.axhline(y = cut_b, color = '#299967', ls = '--')
            ax.set_xlim(med_lims)
            #print "done plot 5"
        elif i == 6:
            ax.plot(np.asarray(med_bw), np.asarray(rms_bw) + np.asarray(med_bw), 'bo')
            ax.set_xlabel('med bw',fontsize=9)
            ax.set_ylabel('rms',fontsize=9)
            ax.set_ylim(rms_lims)
            plt.axhline(y = cut_b, color = '#299967', ls = '--') 
            ax.set_xlim(med_lims)
            #print "done plot 6"    

    plt.savefig(figureName)
    plt.close()
        
        
 
if  __name__ == '__main__':
    main()
