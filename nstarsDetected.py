#!/usr/bin/env python
"""
Extract the number of stars detected for a exposures.
Plot a comparison of number of stars for different nights in order to determine when
do we have a cloudy night.


"""
# imports
import coreutils.desdbi
import cx_Oracle
import argparse
import os
import sys
import numpy as np
import pyfits as pf
import scipy as scp
import os
import math
import re
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.colorbar as cbar
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.mplot3d.axes3d import Axes3D
import time
import collections

parser = argparse.ArgumentParser(description="Extract number of detections per field for all Y2 data.")
parser.add_argument("--section", "-s", required=False, dest="section", help="section of .desservices file with connection info (db-desoper, db-destest)", default="db-desoper")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")


args = parser.parse_args()
section = args.section
verbose = args.verbose

__version__ = "$Revision: 187 $"


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
    cur.execute(query)

    return cur


def main():

    #output filename with summary of data in plot
    summaryOfDataFilename = 'NSTARScut.dat'
    
    field = ['SN-X1','SN-X2','SN-X3','SN-S1','SN-S2','SN-E1','SN-E2','SN-C1','SN-C2','SN-C3']
    bands = ['g','r','i','z']

    dbh = connectDB('db-desoper')
    
    query = """select reqnun,attnum,field from snsubmit where nite >= %s and nite <= %s"""
    
    for myfield in field:
        fig = plt.figure()
        NstarsHSN = {}
        badWeather = {}
        badnites = {}
        #All Exposures That finished Successfully (NEW or Active and Status = 0
        fieldBandNewActive = {}
        #All Exposures in JUNK Runs. JUNK or status != 0
        fieldBandJunk = {}
        #All exposures with HSN < 600
        fieldBandBad = {}

        for myband in bands:
            print myfield,myband
            expnumReqnumAttnum = {}        
        
            mystring = myfield + '_' + myband
            
            #print mystring
            
            query = """SELECT distinct sqa.filename, sqa.astromndets_ref_highsn, sqa.astromndets_ref
                    FROM scamp_qa sqa
                    WHERE sqa.filename LIKE '%%%s%%'
                    order by sqa.filename""" % ( mystring )
                   
            print query
            query_result = query_to_cur(dbh,query)
            
            filename_new = []
            filename_junk = []
            filename_bad = []
            hSN_new = []
            hSN_junk = []
            hSN_bad = []
            filename = []
            ref_new = []
            ref_junk = []
            ref_bad = []
            nite_new = []
            unitname_new = []
            status_new = []
            dataState_new = []
            attnum_new = []
            nite_junk = []
            unitname_junk = []
            status_junk = []
            dataState_junk = []
            attnum_junk = []
            nite_bad = []
            unitname_bad = []
            status_bad = []
            dataState_bad = []
            attnum_bad = []
            the_expnum_new = []
            the_expnum_junk = []
            the_expnum_bad = []
            reqnum_new = []
            reqnum_junk = []
            reqnum_bad = []
            #expnum = []
            #attnum = []
            #reqnum = []
            hSN = []
            ref = []
            
            for val in query_result:
                fileSplitted = val[0].split('_')
                root = fileSplitted[0]
                #if root[1:] in expnum:
                #    continue
                expnum = int(root[1:])
                reqattnum = fileSplitted[3]
                ra = re.match('r(\d+)p(\d+)',reqattnum)
                reqnum = int(ra.group(1))
                attnum = int(ra.group(2))
                hSN = val[1]
                ref = val[2]
                filename = val[0]
                
                #Add to the dictionary only the last attempt for each exposure. Use this to query the database for the rest of the values
                expnumReqnumAttnum[expnum] = [reqnum,attnum,hSN,ref]
                
            
            for expnum,values in expnumReqnumAttnum.items():
                reqnum = values[0]
                attnum = values[1]
                nstarsHSN = values[2]
                ref = values[3]
                
            
            
                #get the nite for each exposure
                #print root
                #get a list of all exposures with sucessful runs (finish with status=0)
                query_snsubmit_new = """select nite,unitname,status,data_state,attnum,reqnum from snsubmit where reqnum=%s and attnum=%s and field='%s' and band='%s' 
                                        and (data_state='NEW' or data_state='ACTIVE') and status=0""" % (reqnum,attnum,myfield,myband)
                #print query_nite
                result_query_snsubmit_new = query_to_cur(dbh, query_snsubmit_new)

                #collect all NEW or active data
                for values in result_query_snsubmit_new:
                    if len(values[0]) > 0:
                        nite_new.append(values[0])
                        unitname_new.append(values[1])
                        status_new.append(values[2])
                        dataState_new.append(values[3])
                        attnum_new.append(values[4])
                        reqnum_new.append(values[5])
                        hSN_new.append(nstarsHSN)
                        ref_new.append(ref)
                        filename_new.append(filename)
                        the_expnum_new.append(expnum)

            
                #query for exposures in JUNK runs and status=0 
                query_snsubmit_junk_success = """select nite,unitname,status,data_state,attnum,reqnum from snsubmit where reqnum=%s and attnum=%s and field='%s' and band='%s'
                                    and data_state = 'JUNK' and status=0""" % (reqnum,attnum,myfield,myband)
                #print query_nite
                result_query_snsubmit_junk_success = query_to_cur(dbh, query_snsubmit_junk_success)

                #Collect all JUNK runs with good status
                for values in result_query_snsubmit_junk_success:
                    if len(values[0]) > 0:
                        nite_junk.append(values[0])
                        unitname_junk.append(values[1])
                        status_junk.append(values[2])
                        dataState_junk.append(values[3])
                        attnum_junk.append(values[4])
                        reqnum_junk.append(values[5])
                        hSN_junk.append(nstarsHSN)
                        ref_junk.append(ref)
                        filename_junk.append(filename)
                        the_expnum_junk.append(expnum)

            
                #query for exposures JUNK and status=1 
                query_snsubmit_junk_fail = """select nite,unitname,status,data_state,attnum,reqnum from snsubmit where reqnum=%s and attnum=%s and field='%s' and band='%s'
                                        and data_state = 'JUNK' and status=1""" % (reqnum,attnum,myfield,myband)
                #print query_nite
                result_query_snsubmit_junk_fail = query_to_cur(dbh, query_snsubmit_junk_fail) 
                
                #Collect all JUNK runs with status=1 (failed)
                for values in result_query_snsubmit_junk_fail:
                    if len(values[0]) > 0:
                        nite_bad.append(values[0])
                        unitname_bad.append(values[1])
                        status_bad.append(values[2])
                        dataState_bad.append(values[3])
                        attnum_bad.append(values[4])
                        reqnum_bad.append(values[5])
                        hSN_bad.append(nstarsHSN)
                        ref_bad.append(ref)
                        filename_bad.append(filename)
                        the_expnum_bad.append(expnum)
                
            
            #All Exposures That finished Successfully (NEW or Active and Status = 0)
            fieldBandNewActive[(myfield,myband)] = [hSN_new,nite_new,unitname_new,attnum_new,status_new,dataState_new,the_expnum_new,reqnum_new]
            #All Exposures in JUNK Runs. JUNK or status = 0
            fieldBandJunk[(myfield,myband)] = [hSN_junk,nite_junk,unitname_junk,attnum_junk,status_junk,dataState_junk,the_expnum_junk,reqnum_junk]
            #All exposures JUNK and status = 1
            fieldBandBad[(myfield,myband)] = [hSN_bad,nite_bad,unitname_bad,attnum_bad,status_bad,dataState_bad,the_expnum_bad,reqnum_bad]

            
            
            #print "keys",fieldBandNewActive.keys()
            #print fieldBandNewActive[(myfield,myband)]
            
        
        
        #print NstarsHSN.keys()
        for i, theband in enumerate(bands):    # i runs from 1 to 6
            xmin = min(fieldBandNewActive[(myfield,theband)][0]+fieldBandJunk[(myfield,theband)][0]+fieldBandBad[(myfield,theband)][0])
            xmax = max(fieldBandNewActive[(myfield,theband)][0]+fieldBandJunk[(myfield,theband)][0]+fieldBandBad[(myfield,theband)][0])
            print xmin,xmax,myfield,theband
            print "NEW :", fieldBandNewActive[(myfield,theband)]
            print "JUNK 0:", fieldBandJunk[(myfield,theband)]
            print "JUNK 1:", fieldBandBad[(myfield,theband)]
            #print "bad weather", badWeather[theband]
            plt.subplot(2,2,i)
            #print theband,i
            
            plt.hist(fieldBandNewActive[(myfield,theband)][0], 30, facecolor='green', alpha=0.9, range=[100,1000])
            legnewact = fieldNewActive[(myfield,theband)][1]
            plt.hist(fieldBandJunk[(myfield,theband)][0],30,facecolor='red', alpha=0.9, range=[xmin,xmax],label=legjunk,linestyle='dotted')
            plt.legend()

            if len(fieldBandJunk[(myfield,theband)][0]) > 0:
                plt.hist(fieldBandJunk[(myfield,theband)][0],30,facecolor='red', alpha=0.6)
            if len(fieldBandBad[(myfield,theband)][0]) > 0:
                plt.hist(fieldBandBad[(myfield,theband)][0],30,facecolor='blue', alpha=0.6)
                
        
            title = myfield + '-' + theband
            plt.title(title)
                
        outfilename = myfield+'.pdf'
        plt.savefig(outfilename)
        plt.close()


        print fieldBandNewActive.keys()
        print fieldBandNewActive.values()
        
        #sort list of each dictionary in expnum order
        
        with open(summaryOfDataFilename,'w') as outfile:
            outfile.write('#FIELD BAND NSTARS EXPNUM NITE DATA_STATE STATUS REQNUM ATTNUM\n')
            #Filename output columns are:
            #field, band, expnum, nite, Data_State, Status, Nstars_HSN, reqnum, attnum 
            #for myfield in field:
            for myband in bands:
                for row in fieldBandNewActive.items():
                    nstars = fieldBandNewActive[(myfield,myband)][0]
                    nites = fieldBandNewActive[(myfield,myband)][1]
                    unitnames = fieldBandNewActive[(myfield,myband)][2]
                    attnums = fieldBandNewActive[(myfield,myband)][3]
                    status = fieldBandNewActive[(myfield,myband)][4]
                    datastates = fieldBandNewActive[(myfield,myband)][5]
                    expnums = fieldBandNewActive[(myfield,myband)][6]
                    reqnums = fieldBandNewActive[(myfield,myband)][7]
                    for i in range(len(nstars)):
                        outfile.write('%s %s %4d %10d %s %10s %2d %5d %2d \n' % (myfield, myband, nstars[i], expnums[i], nites[i], datastates[i], status[i], reqnums[i], attnums[i]))

                for row in fieldBandJunk.items():
                    nstars = fieldBandJunk[(myfield,myband)][0]
                    nites = fieldBandJunk[(myfield,myband)][1]
                    unitnames = fieldBandJunk[(myfield,myband)][2]
                    attnums = fieldBandJunk[(myfield,myband)][3]
                    status = fieldBandJunk[(myfield,myband)][4]
                    datastates = fieldBandJunk[(myfield,myband)][5]
                    expnums = fieldBandJunk[(myfield,myband)][6]
                    reqnums = fieldBandJunk[(myfield,myband)][7]
                    if len(nstars) > 0:
                        for i in range(len(nstars)):
                            outfile.write('%s %s %4d %10d %s %10s %2d %5d %2d \n' % (myfield, myband, nstars[i], expnums[i], nites[i], datastates[i], status[i], reqnums[i], attnums[i]))

                for row in fieldBandBad.items():
                    nstars = fieldBandBad[(myfield,myband)][0]
                    nites = fieldBandBad[(myfield,myband)][1]
                    unitnames = fieldBandBad[(myfield,myband)][2]
                    attnums = fieldBandBad[(myfield,myband)][3]
                    status = fieldBandBad[(myfield,myband)][4]
                    datastates = fieldBandBad[(myfield,myband)][5]
                    expnums = fieldBandBad[(myfield,myband)][6]
                    reqnums = fieldBandBad[(myfield,myband)][7]
                    if len(nstars) > 0:
                        for i in range(len(nstars)):
                            outfile.write('%s %s %4d %10d %s %10s %2d %5d %2d \n' % (myfield, myband, nstars[i], expnums[i], nites[i], datastates[i], status[i], reqnums[i], attnums[i]))

        
                                      
        #plt.show()
    

if  __name__ == '__main__':
    main()
