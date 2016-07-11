#!/usr/bin/env python
"""
For a given night it finds all the observations for SN fields.
It compares the output with a list of "good exposures" given by Marissa
and it creates an .ids file with only the good exposures from Marissa
exposure list file.
It also gets rid off the 10 sec exposures used for pointing.

.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>

the good exposures file given by Marissa is at:

/home/ricardoc/DESDM/Y1A1-SNe/etc/good_expnum_sne_Y1.txt

Example:

createIdList2Process.py -n 20130911

Output will be written at:
/home/ricardoc/DESDM/Y1A1-SNe/snese/ids/

This two directories are hardcoded at this moment to be run for the reprocessing data.


"""

# imports
import coreutils.desdbi
import argparse
import os
import sys
import time
from datetime import datetime, date, time, timedelta
import string
import csv

parser = argparse.ArgumentParser(description="Check the SN Manifest json file sent for every SNe block and when all files arrive submit the firstcut_snse job")
parser.add_argument("-n", "--nite", required=True, dest="night", help="Night to process")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")

args = parser.parse_args()
night = args.night
verbose = args.verbose

    
#    Turn on svn keyword substitution for Revision (above) with:
#    bash$ svn propset svn:keywords "Revision" foo_quality.py
__version__ = "$Revision: 187 $"


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
        print 'which_services_section = ' , dbh.which_services_section()
        
    cur = dbh.cursor()
        
    return cur
               

def create_id_list(night):
    """
    Creates a file with all the exposure.id for all the SNe exposures to process for given night
    :param: night: Night to process
    """
    mycur = connectDB("db-desoper")

    #Open the file with all good exposures
    goodExps = []
    goodExposures = '/home/ricardoc/DESDM/Y1A1-SNe/etc/good_expnum_sne_Y1.txt'
    idsPath = '/home/ricardoc/DESDM/Y1A1-SNe/snese/ids/'

    with open(goodExposures, 'r') as thefile:
        for line in thefile:
            line = line.rstrip()
            exp = line.lstrip()
            exp = exp.rstrip()
            #print "exposure", exp
            goodExps.append(exp)

    
    #print goodExps

    id_query = """select id,expnum,object,exptime,band from exposure
                    where nite = '%s'
                    and object like '%%SN-%%'
                    and project = 'DTS'
                    and exptime > '15'
                    order by exposurename""" % (night)

 
    mycur.arraysize = 1000 # get 1000 at a time when fetching
    mycur.execute(id_query)
    numExps = []
    for item in mycur:
        numExps.append(item)

    totalNumExps = len(numExps)
    print "Total Number of Exposures in this night:", totalNumExps
    print numExps

    shallowFields = ['SN-C1', 'SN-C2', 'SN-X1', 'SN-X2', 'SN-E1', 'SN-E2', 'SN-S1', 'SN-S2']
    deepFields = ['SN-X3', 'SN-C3']
    nextShallowIsGood = False
    nextDeepIsGood = False

    mycur.execute(id_query)
    outfilename = idsPath + 'SN_' + night + '_repro.ids'
    with open(outfilename,'w') as outfile:
        outfile.write('id=')
        for i,item in enumerate(mycur):
            expid = item[0]
            expnum = item[1]
            objectName  = item[2]
            exptime = item[3]
            band = item[4]
            #print expid,objectName,exptime

            for field in shallowFields:
                #print field, objectName
                if field in objectName:
                    isShallow = True
                    isDeep = False
                    thefield = field

            for field in deepFields:
                #print field, objectName
                if field in objectName:
                    isShallow = False
                    isDeep = True
                
            #print field, thefield, isShallow, isDeep, objectName

            if str(expnum) in goodExps:
                if i == totalNumExps - 1:
                    outfile.write(str(expid))
                    print "Good Exposure to process", expnum,exptime,band, objectName
                else:
                    outfile.write(str(expid) + ',')
                    print "Good Exposure to process", expnum,exptime,band, objectName
                if band == 'z' and isShallow:
                    #print "next Exposure is a z band and needs to be processed too"
                    nextShallowIsGood = True
                if isDeep:
                    print 
                    nextDeepIsGood = True
                

            if str(expnum) not in goodExps and isShallow and band == 'z' and nextShallowIsGood:
                #There is one more z band exposure that is not in good list.
                if i == totalNumExps - 1:
                    outfile.write(str(expid))
                    print "Good Exposure to process", expnum,exptime,band, objectName
                else:
                    outfile.write(str(expid) + ',')
                    print "Good Exposure to process", expnum,exptime,band, objectName
                    nextShallowIsGood = False
                continue

            if str(expnum) not in goodExps and isDeep and (band == 'z' or band == 'r' or band == 'i' or band == 'g') and nextDeepIsGood:
                if i == totalNumExps -1:
                    outfile.write(str(expid))
                    print "Good Exposure to process", expnum,exptime,band, objectName
                else:
                    outfile.write(str(expid) + ',')
                    print "Good Exposure to process", expnum,exptime,band, objectName
                continue
            
            if str(expnum) not in goodExps:
                print "The exposure %d - %s - %d - %s is a bad exposure (not present in good exps list)" % (expnum, objectName, exptime, band) 
                
                nextDeepIsGood = False
                nextShallowIsGood = False


#    with open(outfilename, 'r+') as f:
#        allLines = f.readlines()
#        if allLines[-1] == ',':
#            allLines = allLines.rstrip(',')
#            print allLines
        #f.write(allLines)
        
def main():


    create_id_list(night)
    

if  __name__ == '__main__':
    main()


