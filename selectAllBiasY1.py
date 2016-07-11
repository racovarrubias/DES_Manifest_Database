#!/usr/bin/env python
"""Calculates Gain and rdnoise for DES images

.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>

Select all biases taken for year 1, except those present in rgruendl.bias_frame_qa
and those that contain in the object header keyword the word junk or test.
It writes the results in bias2process.txt

  Needs to setup the following before running:
  python
  pyfits
  coreutils
  numpy

"""
# imports
import coreutils.desdbi
import argparse
import os
import sys
import numpy as np
import pyfits as pf
import scipy
import os
import math
import re
import csv


parser = argparse.ArgumentParser(description="Given two raw flats and two raw bias, Calculates gain and rdnoise for all ccds")
parser.add_argument("-n", "--nite", required=False, dest="night", help="Night to select the bias")
parser.add_argument("-t", "--type", required=False, dest="inputData", help="Data type to use for calculation of gain and rdnoise: Options are src or raw")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")

args = parser.parse_args()
night = args.night
verbose = args.verbose
inputData = args.inputData

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
    cur.execute(query)
    #print "Executing query :", query_flats
    #dbcursor.arraysize = 1000 # get 1000 at a time when fetching
    #dbcursor.execute(query_flats)

    return cur



def main():


    #print "band:" ,filter
    #directory where to look for raw data
    if inputData == 'src':
        #src data will be used to calcualte gain and rdnoise
        root_dir = '/archive_data/Archive/DTS/src/'
    elif inputData == 'raw':
        #crostalked, overscan subtracted flats and bias will be used to calculate gain and rdnoise
        root_dir = '/archive_data/Archive/OPS/red/'
    
    dbh = connectDB('db-desoper')


    #queryitems = ["id","expnum","mjd_obs", "nite"]
    queryitems = ["expnum"]
    #querylist = ",".join(queryitems)

    #query_biases = """ select %s from exposure where 
    #                obstype = 'zero' and
    #                project = 'DTS' and
    #                object IS NULL and
    #                band like '%%Empty%%'
    #                order by exposurename""" % (querylist, night)
    
    query_biases = """ select expnum from exposure where 
                    obstype = 'zero' and
                    project = 'DTS' and
                    object NOT LIKE '%%junk%%' and
                    object NOT LIKE '%%test%%' and
                    nite >= '20130815' and nite < '20140220'
                    and expnum not in (select expnum from gruendl.bias_frame_qa)
                    order by expnum"""


    cur = query_to_cur(dbh, query_biases, verbose=True)
    
    
    expnum = []
    mjd_obs = []
    nite = []
    id = []
    for val in cur:
        expnum.append(val[0])
        
    #allnights = set(nite)
    #orderedNites = sorted(allnights)
    
    
    with open("bias2process.txt",'w') as out:
        for exp in sorted(expnum):
            out.write('%d \n' % exp)
    
    sys.exit()
    
    #for each night find all the biases that we should process.
    
    finalExpnums = []
    for singleNight in orderedNites:
        indxSingleNight = [i for i,val in enumerate(nite) if val == singleNight]
        temp_mjdobs = [mjd_obs[i] for i in indxSingleNight]
        temp_expnums = [expnum[i] for i in indxSingleNight]
        temp_expid = [expnum[i] for i in indxSingleNight]
        
        #print "Night = ", singleNight
        #print temp_expnums
        #print temp_mjdobs

    query_bias_frame_qa = """select expnum from exposure where nite > '20130815' and expnum not in 
                            (select expnum from gruendl.bias_frame_qa) order by expnum"""
    
    #query_bias_frame_qa = """select distinct(expnum) from gruendl.bias_frame_qa
    #                        order by expnum"""
                            

    cur_bias = query_to_cur(dbh, query_bias_frame_qa)

    robertTableExpnum = []
    for val in cur_bias:
        robertTableExpnum.append(val[0])
    
    print "Number of bias in bias_frame_qa is", len(robertTableExpnum) 
    print "Number of bias from 20130815", len(expnum)
    print "removing repeated items from robert list", len(set(robertTableExpnum))

    #Eliminate all expnums from all the biases list which are present in robertTable
    #expnumsToProcess = sorted(list(set(robertTableExpnum) - set(expnum)))
    for item in robertTableExpnum:
        if item in expnum:
            expnum.remove(item)
            print "removed exposure from expnum", item
            

    print "new expnum:", len(expnum)
         
    
    #expnumsToProcess = list(set(expnum) - set(robertTableExpnum))
    #allindx = [i for i, val in enumerate(expnum) if val in expnumsToProcess]
    #print allindx
    

    for val in expnum:
        checkExps = """select distinct(%s) from gruendl.bias_frame_qa where """ % (val)
        cur_check = query_to_cur(dbh, checkExps)
        for check in cur_check:
            rexp = check[0]
        if rexp == val:
            print "xxxxx"
        #print val,rexp
        

    sys.exit()
    
    tableExps = []
    for exp in cur_check:
        tableExps.append(exp[0])
    
    
    print len(expnumsToProcess), len(robertTableExpnum), len(expnum)

    
    for val in expnumsToProcess:
        for exp in tableExps:
            print exp, val
            if exp == val:
                print "Exposure is in Robert list ", exp
    
    with open("bias2process.txt",'w') as out:
        for exp in sorted(expnumsToProcess):
            out.write('%d \n' % exp)

    
    #Determines success
    return 0
    
if  __name__ == '__main__':
    status = main()
    sys.exit(status)
    
