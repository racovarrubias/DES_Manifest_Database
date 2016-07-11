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
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
import time
import collections
from pylab import *

#Fontsize for plots
rcParams['axes.labelsize'] = 7
rcParams['xtick.labelsize'] = 7
rcParams['ytick.labelsize'] = 7
rcParams['legend.fontsize'] = 5

#Tick locator
my_locator = MaxNLocator(6)


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

def sortlist(exposures, list2sort):
    """ Sort a the list2sort with repect to the expnum"""
    
    sortedExposures, sortedlist = zip(*sorted(zip(exposures, list2sort)))

    return sortedlist

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    print'axis', x,y,width,height
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    #print "labelsize ", y_labelsize,x_labelsize
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax



def main():

    #output filename with summary of data in plot
    summaryOfDataFilename = 'NSTARScut.dat'
    
    field = ['SN-X1','SN-X2','SN-X3','SN-S1','SN-S2','SN-E1','SN-E2','SN-C1','SN-C2','SN-C3']
    bands = ['g','r','i','z']

    dbh = connectDB('db-desoper')
    
    #query = """select reqnun,attnum,field from snsubmit where nite >= %s and nite <= %s"""

    fieldBandNewActive = {}
    #All Exposures in JUNK Runs. JUNK or status != 0
    fieldBandJunk = {}
    #All exposures with HSN < 600
    fieldBandBad = {}
    fieldBandNewBad = {}
    
    for myfield in field:
        fig = plt.figure()
        NstarsHSN = {}
        badWeather = {}
        badnites = {}
        #All Exposures That finished Successfully (NEW or Active and Status = 0

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
            nite_newbad = []
            unitname_newbad = []
            status_newbad = []
            dataState_newbad = []
            attnum_newbad = []
            reqnum_newbad = []
            hSN_newbad = []
            ref_newbad = []
            filename_newbad = []
            the_expnum_newbad = []

            
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
                #print expnum
            
            
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


                #query for exposures with NEW and status=1 (runs we decided not to process again because of some failure
                #or really bad and cloudy
                query_snsubmit_new_bad = """select nite,unitname,status,data_state,attnum,reqnum from snsubmit where reqnum=%s and attnum=%s and field='%s' and band='%s'
                                    and (data_state = 'NEW' or data_state = 'ACTIVE') and status=1""" % (reqnum,attnum,myfield,myband)
                #print query_nite
                result_query_snsubmit_new_bad = query_to_cur(dbh, query_snsubmit_new_bad)

                #Collect all JUNK runs with good status
                for values in result_query_snsubmit_new_bad:
                    if len(values[0]) > 0:
                        nite_newbad.append(values[0])
                        unitname_newbad.append(values[1])
                        status_newbad.append(values[2])
                        dataState_newbad.append(values[3])
                        attnum_newbad.append(values[4])
                        reqnum_newbad.append(values[5])
                        hSN_newbad.append(nstarsHSN)
                        ref_newbad.append(ref)
                        filename_newbad.append(filename)
                        the_expnum_newbad.append(expnum)

            
                
                
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
                    #print 'yyyyyyyyy', nite_bad, the_expnum_bad
            
            #All Exposures That finished Successfully (NEW or Active and Status = 0)
            fieldBandNewActive[(myfield,myband)] = [hSN_new,nite_new,unitname_new,attnum_new,status_new,dataState_new,the_expnum_new,reqnum_new]
            #All Exposures That failed but were kept because it was not possible to process them (NEW or Active and Status = 1)
            fieldBandNewBad[(myfield,myband)] = [hSN_newbad,nite_newbad,unitname_newbad,attnum_newbad,status_newbad,dataState_newbad,the_expnum_newbad,reqnum_newbad]
            #All Exposures in JUNK Runs. JUNK or status = 0
            fieldBandJunk[(myfield,myband)] = [hSN_junk,nite_junk,unitname_junk,attnum_junk,status_junk,dataState_junk,the_expnum_junk,reqnum_junk]
            #All exposures JUNK and status = 1
            fieldBandBad[(myfield,myband)] = [hSN_bad,nite_bad,unitname_bad,attnum_bad,status_bad,dataState_bad,the_expnum_bad,reqnum_bad]

            
            
            
        
        
    #print NstarsHSN.keys()
        for i, theband in enumerate(bands):    # i runs from 1 to 6
            #xmin = min(fieldBandNewActive[(myfield,theband)][0]+fieldBandJunk[(myfield,theband)][0]+fieldBandBad[(myfield,theband)][0])
            #xmax = max(fieldBandNewActive[(myfield,theband)][0]+fieldBandJunk[(myfield,theband)][0]+fieldBandBad[(myfield,theband)][0])
            #print xmin,xmax,myfield,theband
            #print "NEW :", fieldBandNewActive[(myfield,theband)]
            #print "NEW BAD:, fieldBandNewBad[(myfield,myband)] 
            #print "JUNK 0:", fieldBandJunk[(myfield,theband)]
            #print "JUNK 1:", fieldBandBad[(myfield,theband)]
            #print "bad weather", badWeather[theband]
            ax=plt.subplot(2,2,i)
            print theband,i
            #NEW and status = 0 Green
            if len(fieldBandNewActive[(myfield,theband)][0]) > 0:
               #legnewgood = fieldBandNewActive[(myfield,theband)][1]
                plt.hist(fieldBandNewActive[(myfield,theband)][0], bins=40,color=['green'], alpha=0.9, range=[400,1000])
                title = myfield + '-' + theband
                plt.title(title)
                
            #NEW but failed processing (Yellow)
            if len(fieldBandNewBad[(myfield,theband)][0]) > 0:# and len(fieldBandNewBad[(myfield,theband)][0]) > 0 and len(fieldBandJunk[(myfield,theband)][0]) == 0 and len(fieldBandBad[(myfield,theband)][0]) == 0:
                #xmin = min(fieldBandNewBad[(myfield,theband)][0])
                #xmax = max(fieldBandNewBad[(myfield,theband)][0])
                legnewbad = set(fieldBandNewBad[(myfield,theband)][1])
                print "legnewbad", set(legnewbad)
                selectedNites = []
                selectedIndxs = []
                selectedNstars = []
                otherNites = []
                otherIndxs = []
                otherNstars = []
                #separate data from 20141124 
                if len(legnewbad) > 1:
                    for indx, value in enumerate(fieldBandNewBad[(myfield,theband)][1]):
                        if value == '20141124':
                            selectedNites.append(value)
                            selectedIndxs.append(indx)
                            selectedNstars.append(fieldBandNewBad[(myfield,theband)][0][indx])
                        else:
                            otherNites.append(value)
                            otherIndxs.append(indx)
                            otherNstars.append(fieldBandNewBad[(myfield,theband)][0][indx])               
                    if len(selectedNites) > 0:
                        leg = set(selectedNites)
                        plt.hist(selectedNstars, 40,facecolor='red', alpha=0.6,label=leg,linestyle='dashed', range=[400,1000], rwidth = 2)
                        plt.legend()
                    elif len(otherNites) > 0:
                        leg = set(otherNites)
                        plt.hist(otherNstars, 40,facecolor='yellow', alpha=0.6,label=leg,linestyle='dashed', range=[400,1000], rwidth = 2)
                        plt.legend()
                #else:
                #    plt.hist(fieldBandNewBad[(myfield,theband)][0], 40,facecolor='yellow', alpha=0.6,label=legnewbad,linestyle='dashed', range=[400,1000], rwidth = 2)
                #    plt.legend()
                title = myfield + '-' + theband
                plt.title(title)
                xnum = min(fieldBandNewBad[(myfield,theband)][0])
                #if there is a set of NSTARTS < 400, then make a subplot showing only those values
                if xnum < 400:
                    data_20141124=[]
                    restdataNewBad = []
                    badnite_20141124 = []
                    restBadNiteNewBad = []
                    alldataNewBad = []
                    #axes = plt.figure.add_subplot(2,2,i)
                    xmin = min(fieldBandNewBad[(myfield,theband)][0])
                    position = [0.1, 0.5, .4, .4]
                    subax1 = add_subplot_axes(ax,position)            
                    for indx, val in enumerate(fieldBandNewBad[(myfield,theband)][0]):
                        if val < 400:
                            if fieldBandNewBad[(myfield,theband)][1][indx] == '20141124':
                                data_20141124.append(val)
                                badnite_20141124.append(fieldBandNewBad[(myfield,theband)][1][indx])
                                alldataNewBad.append(fieldBandNewBad[(myfield,theband)][1][indx])
                            else:
                                restdataNewBad.append(value)
                                restBadNiteNewBad.append(fieldBandNewBad[(myfield,theband)][1][indx])
                                alldataNewBad.append(fieldBandNewBad[(myfield,theband)][1][indx])
                                
                    #print myfield, theband, len(x1),x1,badnite
                    max_xticks = 4
                    xloc = plt.MaxNLocator(max_xticks)
                    subax1.xaxis.set_major_locator(xloc)
                    badnites = set(alldataNewBad)
                    for leg in badnites:
                        if leg == '20141124':
                            legend24 = leg                         
                            subax1.hist(data_20141124, 40,facecolor='red', alpha=0.6,linestyle='dashed', range=[xmin-30,400], rwidth = 2,label=legend24)
                            #plt.legend()
                        else:
                            legendrest = set(restBadNiteNewBad)
                            subax1.hist(restdataNewBad, 40,facecolor='yellow', alpha=0.6,linestyle='dashed', range=[xmin-30,400], rwidth = 2,label=legendrest)
                    plt.legend()
                        
                #print myfield, theband, fieldBandNewBad[(myfield,theband)][0]
                
            """
            if len(fieldBandJunk[(myfield,theband)][0]) > 0:
                #xmin = min(fieldBandJunk[(myfield,theband)][0])
                #xmax = min(fieldBandJunk[(myfield,theband)][0])
                legjunk = set(fieldBandJunk[(myfield,theband)][1])
                if len(legjunk) > 1:
                    for leg in legjunk:    
                        plt.hist(fieldBandJunk[(myfield,theband)][0],40,facecolor='red', alpha=0.6,label=leg,linestyle='dotted', range=[400,1000],rwidth = 2)
                        plt.legend()
                else:
                    plt.hist(fieldBandJunk[(myfield,theband)][0],40,facecolor='red', alpha=0.6,label=legjunk,linestyle='dotted', range=[400,1000],rwidth = 2)
                    plt.legend()
                title = myfield + '-' + theband
                plt.title(title)
                xnum = min(fieldBandJunk[(myfield,theband)][0])
                #if there is a set of NSTARTS < 400, then make a subplot showing only those values
                if xnum < 400:
                    x2=[]
                    #axes = plt.figure.add_subplot(2,2,i)
                    position = [0.1, 0.5, .4, .4]
                    xmin = min(fieldBandJunk[(myfield,theband)][0])
                    subax1 = add_subplot_axes(ax,position)
                    for val in fieldBandJunk[(myfield,theband)][0]:
                        if val < 400:
                            x2.append(val)
                        print len(x2),x2
                    print myfield, theband, len(x2),x2
                    max_xticks = 4
                    xloc = plt.MaxNLocator(max_xticks)
                    subax1.xaxis.set_major_locator(xloc)

                    subax1.hist(x2, 40,facecolor='red', alpha=0.6,linestyle='dashed', range=[xmin-30,400], rwidth = 2)
            """
            #JUNK and failed processing staus = 1 --> Blue
            if len(fieldBandBad[(myfield,theband)][0]) > 0:
            ##    xmin = min(fieldBandBad[(myfield,theband)][0])
            #   xmax = min(fieldBandBad[(myfield,theband)][0])
                legbad = set(fieldBandBad[(myfield,theband)][1])
                selectedNites = []
                selectedIndxs = []
                selectedNstars = []
                otherNites = []
                otherIndxs = []
                otherNstars = []
                #separate data from 20141124 
                if len(legbad) > 1:
                    for indx, value in enumerate(fieldBandBad[(myfield,theband)][1]):
                        if value == '20141124':
                            selectedNites.append(value)
                            selectedIndxs.append(indx)
                            selectedNstars.append(fieldBandBad[(myfield,theband)][0][indx])
                        else:
                            otherNites.append(value)
                            otherIndxs.append(indx)
                            otherNstars.append(fieldBandBad[(myfield,theband)][0][indx])               
                    if len(selectedNites) > 0:
                        leg = set(selectedNites)
                        plt.hist(selectedNstars, 40,facecolor='red', alpha=0.6,label=leg,linestyle='dashed', range=[400,1000], rwidth = 2)
                        plt.legend()
                    elif len(otherNites) > 0:
                        leg = set(otherNites)
                        plt.hist(otherNstars, 40,facecolor='yellow', alpha=0.6,label=leg,linestyle='dashed', range=[400,1000], rwidth = 2)
                        plt.legend()

                title = myfield + '-' + theband
                plt.title(title)
                xnum = min(fieldBandBad[(myfield,theband)][0])
                #if there is a set of NSTARTS < 400, then make a subplot showing only those values
                if xnum < 400:
                    badnites = []
                    data_20141124=[]
                    restdata = []
                    badnite_20141124 = []
                    restBadNite = []
                    alldataBadBad = []
                    #axes = plt.figure.add_subplot(2,2,i)
                    postion = [0.1, 0.5, .4, .4]
                    xmin = min(fieldBandBad[(myfield,theband)][0])
                    subax1 = add_subplot_axes(ax,position)
                    for indx, val in enumerate(fieldBandBad[(myfield,theband)][0]):
                        if val < 400:
                            if fieldBandBad[(myfield,theband)][1][indx] == '20141124':
                                data_20141124.append(val)
                                badnite_20141124.append(fieldBandBad[(myfield,theband)][1][indx])
                                alldataBadBad.append(fieldBandBad[(myfield,theband)][1][indx])
                            else:
                                restdata.append(value)
                                restBadNite.append(fieldBandBad[(myfield,theband)][1][indx])
                                alldataBadBad.append(fieldBandBad[(myfield,theband)][1][indx])
                                
                    #print myfield, theband, len(x3),x3,badnites
                    max_xticks = 4
                    xloc = plt.MaxNLocator(max_xticks)
                    subax1.xaxis.set_major_locator(xloc)
                    legnites = set(alldataBadBad)
                    for leg in legnites:
                        if leg == '20141124':
                            legend24 = leg
                            subax1.hist(data_20141124, 40,facecolor='red', alpha=0.8,linestyle='dashed', range=[xmin-30,400], rwidth = 2, label=legend24)
                        else:
                            legendrest = set(restBadNite)
                            subax1.hist(restdata, 40,facecolor='blue', alpha=0.8,linestyle='dashed', range=[xmin-30,400], rwidth = 2, label=legendrest)
                    plt.legend()
                            
        outfilename = myfield+'.pdf'
        plt.savefig(outfilename)
        plt.close()


        #print fieldBandNewActive.keys()
        #print fieldBandNewActive.values()
        
        #sort list of each dictionary in expnum order
        
        with open(summaryOfDataFilename,'a') as outfile:
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
                    #sort all list by expnum
                    nstarsSorted = sortlist(expnums, nstars)
                    nitesSorted = sortlist(expnums, nites)
                    unitnamesSorted = sortlist(expnums, unitnames)
                    attnumsSorted = sortlist(expnums, attnums)
                    statusSorted = sortlist(expnums, status)
                    datastatesSorted = sortlist(expnums, datastates)
                    reqnumsSorted = sortlist(expnums, reqnums)
                    expnumsSorted = sorted(expnums)
        
                    for i in range(len(nstars)):
                        outfile.write('%s %s %4d %10d %s %10s %2d %5d %2d \n' % (myfield, myband, nstarsSorted[i], expnumsSorted[i], nitesSorted[i], datastatesSorted[i], statusSorted[i], reqnumsSorted[i], attnumsSorted[i]))
        
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
                        nstarsSorted = sortlist(expnums, nstars)
                        nitesSorted = sortlist(expnums, nites)
                        unitnamesSorted = sortlist(expnums, unitnames)
                        attnumsSorted = sortlist(expnums, attnums)
                        statusSorted = sortlist(expnums, status)
                        datastatesSorted = sortlist(expnums, datastates)
                        reqnumsSorted = sortlist(expnums, reqnums)
                        expnumsSorted = sorted(expnums)
                        for i in range(len(nstars)):
                            outfile.write('%s %s %4d %10d %s %10s %2d %5d %2d \n' % (myfield, myband, nstarsSorted[i], expnumsSorted[i], nitesSorted[i], datastatesSorted[i], statusSorted[i], reqnumsSorted[i], attnumsSorted[i]))
        
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
                        nstarsSorted = sortlist(expnums, nstars)
                        nitesSorted = sortlist(expnums, nites)
                        unitnamesSorted = sortlist(expnums, unitnames)
                        attnumsSorted = sortlist(expnums, attnums)
                        statusSorted = sortlist(expnums, status)
                        datastatesSorted = sortlist(expnums, datastates)
                        reqnumsSorted = sortlist(expnums, reqnums)
                        expnumsSorted = sorted(expnums)
                        for i in range(len(nstars)):
                            outfile.write('%s %s %4d %10d %s %10s %2d %5d %2d \n' % (myfield, myband, nstarsSorted[i], expnumsSorted[i], nitesSorted[i], datastatesSorted[i], statusSorted[i], reqnumsSorted[i], attnumsSorted[i]))

        
                                      
        #plt.show()
    

if  __name__ == '__main__':
    main()
