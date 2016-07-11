#!/usr/bin/env python

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
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.colorbar as cbar
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.mplot3d.axes3d import Axes3D


parser = argparse.ArgumentParser(description="Check if all flats from a night were processed a")
parser.add_argument("-f", "--file", required=False, dest="inputfile", help="Filename with nites to check")
parser.add_argument("--section", '-s', required=False, dest="section", help='section of .desservices file with connection info (db-desoper, db-destest)', default="db-desoper")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")



args = parser.parse_args()
inputfile = args.inputfile
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
    
    

    dbh = connectDB(section)
    
    with open('nights_to_process.txt', 'w') as outfile:                
        with open(inputfile,'r') as infile:
            for line in infile:
                nite = line.rstrip()
                #print nite
        
                checkNite = """select count(id) from exposure where obstype='dome flat' and nite = '%s' """ % nite
                cursorIn = query_to_cur(dbh,checkNite,section)
                checkIngested = """select count(distinct(exposureid)) from flats_stats_qa where nite= '%s' """ % nite
                cursorIngested = query_to_cur(dbh,checkIngested, section)
                
                for vals in cursorIn:
                    print "value input", vals
                    input = vals[0]
                for val in cursorIngested:
                    print "value ingested", val
                    ingested = val[0]
                    
                print input, ingested, nite    
                if input - ingested != 0:
                    print "NIGHT NOT FULLY PROCESSED: ", nite
                    out = '%s %d %d \n' % (nite, input, ingested)
                    outfile.write(out)
    
    return 0

if  __name__ == '__main__':
    status = main()
    sys.exit(status)
    