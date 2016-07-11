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


parser = argparse.ArgumentParser(description="Examine bias from gruendl.bias_frames_qa_v2 table and assess them")
parser.add_argument("-db", "--insertdb", required=False,action="store_true",default=False, help="Insert results into the database (Default is False). If option, then True")
parser.add_argument("-v", "--verbose", required=False,action="store_true",default=False,help="Use --verbose for verbose output")

args = parser.parse_args()
insertIntoDb = args.insertdb
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

def medsigclip(indata, clipsig=3.0, maxiter=10, converge_num=0.001, verbose=verbose):
   """
   Computes an iteratively sigma-clipped median on a
   data set. Clipping is done about median, but media and mean
   are returned.
  
   :History:
       * 21/10/1998 Written by RSH, RITSS
       * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
       * 24/11/2009 Converted to Python. PLL.
 
   Examples
   --------
   >>> median, mean, sigma = meanclip(indata)
 
   Parameters
   ----------
   indata: array_like
       Input data.
 
   clipsig: float
       Number of sigma at which to clip.
 
   maxiter: int
       Ceiling on number of clipping iterations.
 
   converge_num: float
       If the proportion of rejected pixels is less than
       this fraction, the iterations stop.
 
   verbose: {0, 1}
       Print messages to screen?
 
   Returns
   -------
   median: float
       N-sigma clipped median.

   mean: float
       N-sigma clipped mean.
 
   sigma: float
       Standard deviation of remaining pixels.
 
   """
   # Flatten array
   skpix = indata.reshape( indata.size, )
 
   ct = indata.size
   iter = 0; c1 = 1.0 ; c2 = 0.0
 
   while (c1 >= c2) and (iter < maxiter):
       lastct = ct
       medval = np.median(skpix)
       sig = np.std(skpix)
       wsm = np.where( abs(skpix-medval) < clipsig*sig )
       ct = len(wsm[0])
       if ct > 0:
           skpix = skpix[wsm]
 
       c1 = abs(ct - lastct)
       c2 = converge_num * lastct
       iter += 1
   # End of while loop
 
   mean  = np.mean( skpix )
   sigma = np.std( skpix )
   median = np.median( skpix )
   
   if verbose:
       prf = 'MEDSIGCLIP:'
       print '%s %.1f-sigma clipped mean' % (prf, clipsig)
       print '%s Mean computed in %i iterations' % (prf, iter)
       print '%s , Median = %.6f, Mean = %.6f, sigma = %.6f' % (prf, medval, mean, sigma)
 
   return median, mean, sigma


def insert2Db(dbh, query, verbose=verbose):
    """Execute a query and return a cursor to a query
    :param dbh: Database connection handler
    :param query: query to execute
    :param debug: verbosity
    
    """
    
    print "inside"
    try:      
        #dbh = connectDB('db-desoper')
        cur = dbh.cursor()
        cur.execute(query)
        dbh.commit()
        if verbose:
            print query
    except cx_Oracle.IntegrityError as e:
        print "error while inserting into FLATS_STATS_QA TABLE: ", e  

    #dbh.commit()

def main():


    """
    gruendl.bias_frame_qa
     Name                            Null?    Type
     ----------------------------------------- -------- ----------------------------
     EXPNUM                              NUMBER(10)
     CCD                                  NUMBER(2)
     SOURCE                        NOT NULL VARCHAR2(30)
     MED_A                                  NUMBER(10,3)
     MED_AR                              NUMBER(10,3)
     MED_AW                              NUMBER(10,3)
     MED_B                                  NUMBER(10,3)
     MED_BR                              NUMBER(10,3)
     MED_BW                              NUMBER(10,3)
     RMS_A                                  NUMBER(10,3)
     RMS_AR                              NUMBER(10,3)
     RMS_AW                              NUMBER(10,3)
     RMS_B                                  NUMBER(10,3)
     RMS_BR                              NUMBER(10,3)
     RMS_BW                              NUMBER(10,3)
    
    """

    section = 'db-desoper'
    dbh = connectDB(section)

    allccds = range(1,61)
    allccds.append(62)
    
    statsBias = {}
    allStatsBias = []
    badBias = {}
    goodBias = {}
    
    bkp1 = [4,8,13,14,19,20,25,26,27,32,33,34,39,40,45,46,51,56]
    bkp3 = [1,2,3,5,6,7,9,10,11,12,15,16,17,18,21,22,23,24]
    bkp4 = [28,29,30,31,35,36,37,38]
    bkp5 = [41,42,43,44,47,48,49,50,52,53,54,55,57,58,59,60,61,62]
    
    
    bad_rms_a = []
    bad_med_a = []
    good_rms_a = []
    godd_med_a = []

    bad_rms_b = []
    bad_med_b = []
    good_rms_b = []
    godd_med_b = []

    

    for ccdnum in allccds:

        select_data = """select expnum, ccd, med_a, 
                        med_ar,med_aw, rms_a, rms_ar, 
                        rms_aw, med_b, med_br,med_bw, 
                        rms_b, rms_br, rms_bw, source 
                        from gruendl.bias_frame_qa_v2 
                        where ccd = %s  and 
                        expnum not in (select expnum 
                        from exposure where obstype = 'zero' and
                        object like '%%PTC%%' and
                        object like '%%VSUB%%') """ % (ccdnum)
    
    
        data = query_to_cur(dbh, select_data, verbose=False)
    
        #print rsq_chi2        
        ccd = []
        expnum = []
        med_a = []
        med_ar = []
        med_aw = []
        rms_a = []
        rms_ar = []
        rms_aw = []
        med_b = []
        med_br = []
        med_bw = []
        rms_b = []
        rms_br = []
        rms_bw = []
        source = []
        
        for val in data:
            expnum.append(val[0])
            ccd.append(val[1])
            med_a.append(val[2])
            med_ar.append(val[3])
            med_aw.append(val[4])
            rms_a.append(val[5])
            rms_ar.append(val[6])
            rms_aw.append(val[7])
            med_b.append(val[8])
            med_br.append(val[9])
            med_bw.append(val[10])
            rms_b.append(val[11])
            rms_br.append(val[12])
            rms_bw.append(val[13])
            source.append(val[14])
        
        
        #print ccdnum, len(med_a), len(med_b), len(rms_a), len(rms_b)
        #calculate the mean of both amps outsind sigma clipping rejection
        median_A, mean_A, sigma_A = medsigclip(np.asarray(med_a), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)
        median_AR, mean_AR, sigma_AR = medsigclip(np.asarray(med_ar), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)
        median_AW, mean_AW, sigma_AW = medsigclip(np.asarray(med_aw), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)
        median_B, mean_B, sigma_B = medsigclip(np.asarray(med_b), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)
        median_BR, mean_BR, sigma_BR = medsigclip(np.asarray(med_br), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)
        median_BW, mean_BW, sigma_BW = medsigclip(np.asarray(med_bw), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)

        medianRMS_A, meanRMS_A, sigmaRMS_A = medsigclip(np.asarray(rms_a), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)
        medianRMS_AR, meanRMS_AR, sigmaRMS_AR = medsigclip(np.asarray(rms_ar), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)
        medianRMS_AW, meanRMS_AW, sigmaRMS_AW = medsigclip(np.asarray(rms_aw), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)
        medianRMS_B, meanRMS_B, sigmaRMS_B = medsigclip(np.asarray(rms_b), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)
        medianRMS_BR, meanRMS_BR, sigmaRMS_BR = medsigclip(np.asarray(rms_br), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)
        medianRMS_BW, meanRMS_BW, sigmaRMS_BW = medsigclip(np.asarray(rms_bw), clipsig=2, maxiter=10, converge_num=0.001, verbose=False)

        statsBias['CCDNUM'] = ccdnum
        statsBias['median_A'] = median_A
        statsBias['mean_A'] = mean_A
        statsBias['sigma_A'] = sigma_A
        statsBias['median_AR'] = median_AR
        statsBias['mean_AR'] = mean_AR
        statsBias['sigma_AR'] = sigma_AR
        statsBias['median_AW'] = median_AW
        statsBias['mean_AW'] = mean_AW
        statsBias['sigma_AW'] = sigma_AW
        
        statsBias['median_B'] = median_B
        statsBias['mean_B'] = mean_B
        statsBias['sigma_B'] = sigma_B
        statsBias['median_BR'] = median_BR
        statsBias['mean_BR'] = mean_BR
        statsBias['sigma_BR'] = sigma_BR
        statsBias['median_BW'] = median_BW
        statsBias['mean_BW'] = mean_BW
        statsBias['sigma_BW'] = sigma_BW


        #if ccdnum == 45:
        #    n, bins, patches = plt.hist(rms_b, 50, normed=1, histtype='stepfilled')
            #P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        #    plt.show()
            
        #List of dictionaries with all the values calculated.
        allStatsBias.append(statsBias)

        if verbose:
            
            print "CCD = %2d : Median A  %10.2f, Mean A  %10.2f, Sigma A  %10.2f" % (ccdnum, median_A, mean_A, sigma_A)
            print "CCD = %2d : Median AR %10.2f, Mean AR %10.2f, Sigma AR %10.2f" % (ccdnum, median_AR, mean_AR, sigma_AR)
            print "CCD = %2d : Median AW %10.2f, Mean AW %10.2f, Sigma AW %10.2f" % (ccdnum, median_AW, mean_AW, sigma_AW)
            print "CCD = %2d : Median B  %10.2f, Mean B  %10.2f, Sigma B  %10.2f" % (ccdnum, median_B, mean_B, sigma_B)
            print "CCD = %2d : Median BR %10.2f, Mean BR %10.2f, Sigma BR %10.2f" % (ccdnum, median_BR, mean_BR, sigma_BR)
            print "CCD = %2d : Median BW %10.2f, Mean BW %10.2f, Sigma BW %10.2f" % (ccdnum, median_BW, mean_BW, sigma_BW)
    
            print "CCD = %2d : Median RMS_A  %10.2f, Mean RMS_A  %10.2f, Sigma RMS_A  %10.2f" % (ccdnum, medianRMS_A, meanRMS_A, sigmaRMS_A)
            print "CCD = %2d : Median RMS_AR %10.2f, Mean RMS_AR %10.2f, Sigma RMS_AR %10.2f" % (ccdnum, medianRMS_AR, meanRMS_AR, sigmaRMS_AR)
            print "CCD = %2d : Median RMS_AW %10.2f, Mean RMS_AW %10.2f, Sigma RMS_AW %10.2f" % (ccdnum, medianRMS_AW, meanRMS_AW, sigmaRMS_AW)
            print "CCD = %2d : Median RMS_B  %10.2f, Mean RMS_B  %10.2f, Sigma RMS_B  %10.2f" % (ccdnum, medianRMS_B, meanRMS_B, sigmaRMS_B)
            print "CCD = %2d : Median RMS_BR %10.2f, Mean RMS_BR %10.2f, Sigma RMS_BR %10.2f" % (ccdnum, medianRMS_BR, meanRMS_BR, sigmaRMS_BR)
            print "CCD = %2d : Median RMS_BW %10.2f, Mean RMS_BW %10.2f, Sigma RMS_BW %10.2f" % (ccdnum, medianRMS_BW, meanRMS_BW, sigmaRMS_BW)


        
        #cut_a = median_A + medianRMS_A
        #cut_b = median_B + medianRMS_B
        cut_a = median_A + 2*medianRMS_A
        cut_b = median_B + 2*medianRMS_B

        rms_cut_a = medianRMS_A + medianRMS_A
        rms_cut_b = medianRMS_B + medianRMS_B

            
        print "CCDNUM = %s, cut_A = %3.1f, cut_B = %3.1f" % (ccdnum, cut_a, cut_b)
        indx_a = [i for i,val in enumerate(med_a) if val+rms_a[i] > cut_a]
        indx_b = [i for i,val in enumerate(med_b) if val+rms_b[i] > cut_b]

        #make a plot for each ccd        
        plt.clf()
        #plt.xticks(rotation=45)
        plt.grid()
        plt.rc("font", size=9)

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

        if ccdnum in bkp1:
            figureName = 'Median_AB_Bias_bkp1_' + str(ccdnum) + '.png'
        if ccdnum in bkp3:
            figureName = 'Median_AB_Bias_bkp3_' + str(ccdnum) + '.png'
        if ccdnum in bkp4:
            figureName = 'Median_AB_Bias_bkp4_' + str(ccdnum) + '.png'
        if ccdnum in bkp5:
            figureName = 'Median_AB_Bias_bkp5_' + str(ccdnum) + '.png'
        plt.savefig(figureName)
        plt.close()
        
        
        #The cuts I am applying for each back plane are:
        #Bkp1:
        #3.3 or 4 (all these is about 2*Median RMS
        #Bkp3:
        #3.3
        #Bkp4:
        #Try 3.3 and see what I get. Most are around 4
        #Bkp5:
        #3.3 tb.
        
        
    
    
        #select all bias that have the following cuts:
        #if rms > 3.3
        #Look for all indexes and find the expnum and ccd where rms > 3.3
        #bad_indx_a = [indx for indx, val in enumerate(rms_a) if val >= 3.3] 
        #bad_indx_b = [indx for indx, val in enumerate(rms_b) if val >= 3.3] 
        #good_indx_a = [indx for indx, val in enumerate(rms_a) if val <= 3.3] 
        #good_indx_b = [indx for indx, val in enumerate(rms_b) if val <= 3.3] 


        bad_flag_expnum_A = []
        bad_flag_expnum_B = []
        bad_only_expnum_A = []
        bad_only_expnum_B = []
        bad_flag_ccds_A = []
        bad_flag_ccds_B = []
        bad_flag_source_A = []
        bad_flag_source_B = []
    
        good_flag_expnum_A = []
        good_flag_expnum_B = []
        good_only_expnum_A = []
        good_only_expnum_B = []
        good_flag_ccds_A = []
        good_flag_ccds_B = []
        good_flag_source_A = []
        good_flag_source_B = []



        bad_indx_a = [indx for indx, val in enumerate(med_a) if val+rms_a[indx] >= cut_a]
        bad_indx_b = [indx for indx, val in enumerate(med_b) if val+rms_b[indx] >= cut_b] 
        good_indx_a = [indx for indx, val in enumerate(med_a) if (val+rms_a[indx] < cut_a and val > 0)] 
        good_indx_b = [indx for indx, val in enumerate(med_b) if (val+rms_b[indx] < cut_b and val > 0)] 


        #Print all expnum,ccd,band for biases that have a mean < 0
        indx_neg_a = [indx for indx, val in enumerate(med_a) if val+rms_a[indx] <= 0.0 ]
        indx_neg_b = [indx for indx, val in enumerate(med_b) if val+rms_b[indx] <= 0.0 ]
        
        
        """
        for i in bad_indx_a:
            print "BAD A: expnum %d med %f and val+rms = %f" % (expnum[i], med_a[i], med_a[i]+rms_a[i])
        
        for i in bad_indx_b:
            print "BAD B: expnum %d med %f and val+rms = %f" % (expnum[i], med_b[i], med_b[i]+rms_b[i])

        for i in bad_indx_a:
            print "NEG A: expnum %d med %f and val+rms = %f" % (expnum[i], med_a[i], med_a[i]+rms_a[i])
        
        for i in bad_indx_b:
            print "NEG B: expnum %d med %f and val+rms = %f" % (expnum[i], med_b[i], med_b[i]+rms_b[i])


        for i in good_indx_a:
            print "good A: expnum %d med %f and val+rms = %f" % (expnum[i], med_a[i], med_a[i]+rms_a[i])
        
        for i in good_indx_b:
            print "good B: expnum %d med %f and val+rms = %f" % (expnum[i], med_b[i], med_b[i]+rms_b[i])

        for i,med in enumerate(med_a):
            print "ALL A: expnum %d med %d and val+rms = %f" % (expnum[i], med, med+rms_a[i])
        """
        
        impath = '/archive_data/Archive/OPS/red/'
        decam = 'DECam_00'
     
        #print len(expnum), len(ccd), len(med_a)
        for i in bad_indx_a:        
            #print "A:", expnum[i], ccd[i]
            bad_only_expnum_A.append(expnum[i])
            full_expnum = impath + source[i] +'/raw/' + decam + str(expnum[i]) + '/' + decam + str(expnum[i]) +'_' + str(ccd[i]) + '.fits.fz'
            bad_flag_expnum_A.append(full_expnum)
            bad_flag_ccds_A.append(ccd[i])
            bad_flag_source_A.append(source[i])
            bad_rms_a.append(rms_a[i])
            bad_med_a.append(med_a[i])
            #insert_statement = """insert into ricardoc.bias_qa (id,expnum, ccd, bias_status_a) values (bias_qa_seq.nextval, %d, %s, VETO)""" % (expnum[i], ccdnum)
            #print "STATETEM", insert_statement
            
        for i in bad_indx_b:        
            #print 'B:', expnum[i], ccd[i]
            bad_only_expnum_B.append(expnum[i])
            full_expnum = impath + source[i] +'/raw/' + decam + str(expnum[i]) + '/' + decam + str(expnum[i]) +'_' + str(ccd[i]) + '.fits.fz'
            bad_flag_expnum_B.append(full_expnum)
            bad_flag_ccds_B.append(ccd[i])
            bad_flag_source_B.append(source[i])
            bad_rms_b.append(rms_b[i])
            bad_med_b.append(med_b[i])
            #insert_statement = """insert into ricardoc.bias_qa (bias_status_b) values (%s) where ccd=%s and expnum=%s""" % ('VETO', ccdnum, expnum[i])

        for i in indx_neg_a:        
            #print "A:", expnum[i], ccd[i]
            bad_only_expnum_A.append(expnum[i])
            full_expnum = impath + source[i] +'/raw/' + decam + str(expnum[i]) + '/' + decam + str(expnum[i]) +'_' + str(ccd[i]) + '.fits.fz'
            bad_flag_expnum_A.append(full_expnum)
            bad_flag_ccds_A.append(ccd[i])
            bad_flag_source_A.append(source[i])
            bad_rms_a.append(rms_a[i])
            bad_med_a.append(med_a[i])

        for i in indx_neg_b:        
            #print 'B:', expnum[i], ccd[i]
            bad_only_expnum_B.append(expnum[i])
            full_expnum = impath + source[i] +'/raw/' + decam + str(expnum[i]) + '/' + decam + str(expnum[i]) +'_' + str(ccd[i]) + '.fits.fz'
            bad_flag_expnum_B.append(full_expnum)
            bad_flag_ccds_B.append(ccd[i])
            bad_flag_source_B.append(source[i])
            bad_rms_b.append(rms_b[i])
            bad_med_b.append(med_b[i])


        for i in good_indx_a:        
            #print "A:", expnum[i], ccd[i]
            good_only_expnum_A.append(expnum[i])
            full_expnum = impath + source[i] +'/raw/' + decam + str(expnum[i]) + '/' + decam + str(expnum[i]) +'_' + str(ccd[i]) + '.fits.fz'
            good_flag_expnum_A.append(full_expnum)
            good_flag_ccds_A.append(ccd[i])
            good_flag_source_A.append(source[i])
            good_rms_a.append(rms_a[i])
            godd_med_a.append(med_a[i])


        for i in good_indx_b:        
            #print 'B:', expnum[i], ccd[i]
            good_only_expnum_B.append(expnum[i])
            full_expnum_B = impath + source[i] +'/raw/' + decam + str(expnum[i]) + '/' + decam + str(expnum[i]) +'_' + str(ccd[i]) + '.fits.fz'
            good_flag_expnum_B.append(full_expnum)
            good_flag_ccds_B.append(ccd[i])
            good_flag_source_B.append(source[i])
            good_rms_a.append(rms_a[i])
            godd_med_a.append(med_a[i])



        print "CCDNUM = %s, GOOD_A = %d, GOODB_B = %d" % (ccdnum, len(good_indx_a), len(good_indx_b))
        print "CCDNUM = %s, BAD_A = %d, BAD_B = %d" % (ccdnum, len(bad_indx_a), len(bad_indx_b))
        print "CCDNUM = %s, ALL_A = %d, ALL_B = %d" % (ccdnum, len(med_a), len(med_b))
        print "CCDNUM = %s, NEG_A = %d, NEG_B = %d" % (ccdnum, len(indx_neg_a), len(indx_neg_b))
        
        
        #Populate good and bad bias dictionaries
        
        if insertIntoDb:     
            #merge into ricardoc.exposure_qa exqa using (select id from exposure) eid on (eid.id = exqa.exposureid) when matched then update set exqa.MSURTEMP = 40.0;
            for i,exp in enumerate(bad_only_expnum_A):
                #insert_statement_A = """insert into ricardoc.bias_qa_v2 (id,expnum, ccd, bias_status_a) values (bias_qa_v2_seq.nextval,%d, %d, 'VETO')""" % (exp, bad_flag_ccds_A[i])  
                update_statement_A = """update ricardoc.bias_qa_v2 SET bias_status_a = 'VETO' where expnum = %d and ccd = %d"""% (exp, bad_flag_ccds_A[i])  
                #print "YYYY", insert_statement_A
                insert2Db(dbh, update_statement_A,verbose=False)
            for i,exp in enumerate(good_only_expnum_A):
                #insert_statement_A = """insert into ricardoc.bias_qa_v2 (id,expnum, ccd, bias_status_a) values (bias_qa_v2_seq.nextval,%d, %d, 'GOOD')""" % (exp, good_flag_ccds_A[i])  
                update_statement_A = """update ricardoc.bias_qa_v2 SET bias_status_a = 'GOOD' where expnum = %d and ccd = %d"""% (exp, good_flag_ccds_A[i])  
                insert2Db(dbh, update_statement_A,verbose=False)
                #print "YYYY", insert_statement_A
            for i,exp in enumerate(bad_only_expnum_B):
                insert_statement_B = """update ricardoc.bias_qa_v2 SET bias_status_b = 'VETO' where expnum = %d and ccd = %d""" % (exp, bad_flag_ccds_B[i]) 
                insert2Db(dbh, insert_statement_B,verbose=False)
                #print "XXXX", insert_statement_B
            for i,exp in enumerate(good_only_expnum_B):
                insert_statement_B = """update ricardoc.bias_qa_v2 SET bias_status_b = 'GOOD' where expnum = %d and ccd = %d""" % (exp, good_flag_ccds_B[i]) 
                insert2Db(dbh, insert_statement_B,verbose=False)
                #print "XXXX", insert_statement_B
        
    
    #create a table with veto and good bias: Primary key is expnum
    #ricardoc.bias_good_veto
    #if table exists, then delete it and create it again
    delete_table_query = """drop table ricardoc.bias_good_veto purge"""
    insert2Db(dbh, delete_table_query,verbose=True)    

    #create table
    create_table = """create table ricardoc.bias_good_veto (
                    expnum NUMBER(10) NOT NULL,
                    BIAS_STATUS VARCHAR2(25),
                    constraint bias_good_veto_pk PRIMARY KEY (expnum)) """
                    
    insert2Db(dbh, create_table, verbose=True)
    
    

    #Select all expnum that has at least one VETO status_bias_a or b and declare it as VETO completely.
    #count how many ccd are VETO per expnum
    
    if not insertIntoDb:
        #Select all the exposures that contains VETO. and write them into a file
        #select b.expnum,b.ccd,b.bias_status_a,b.bias_status_b,e.nite from bias_qa_v2 b, exposure e  where b.bias_status_a = 'VETO' 
        #and b.bias_status_b = 'VETO' and e.expnum = b.expnum and e.nite > 20130907 and e.nite < 20140207  order by b.expnum,b.ccd;
        vetoVetoQuery = """select b.expnum,b.ccd,b.bias_status_a,b.bias_status_b,e.nite from bias_qa_v2 b, exposure e  where b.bias_status_a = 'VETO' 
                    and b.bias_status_b = 'VETO' and e.expnum = b.expnum and e.nite > 20130801 and e.nite < 201402010  and b.ccd !=2 order by b.expnum,b.ccd """
        
        vetoGoodQuery = """select b.expnum,b.ccd,b.bias_status_a,b.bias_status_b,e.nite from bias_qa_v2 b, exposure e  where 
                    (b.bias_status_a = 'VETO' and b.bias_status_b = 'GOOD') 
                    and e.expnum = b.expnum and e.nite > 20130801 and e.nite < 201402010  and b.ccd != 2 order by b.expnum,b.ccd """
        
        goodVetoQuery = """select b.expnum,b.ccd,b.bias_status_a,b.bias_status_b,e.nite from bias_qa_v2 b, exposure e  where 
                    (b.bias_status_a = 'GOOD' and b.bias_status_b = 'VETO')  
                    and e.expnum = b.expnum and e.nite > 20130801 and e.nite < 201402010  and b.ccd != 2 order by b.expnum,b.ccd """

        
        vetoVeto = query_to_cur(dbh, vetoVetoQuery, verbose=verbose)
        vetoGood = query_to_cur(dbh, vetoGoodQuery, verbose=verbose)
        goodVeto = query_to_cur(dbh, goodVetoQuery, verbose=verbose)
        
        vetoVetoExpnum = []
        vetoVetoCcd = []
        statvetoVetoA = []
        statvetoVetoB = []
        niteVetoVeto = []
        with open('vetoveto.txt', 'w') as veto:
            for vals in vetoVeto:
                veto.write("%10d %2d %10s %7s %7s \n" % (vals[0],vals[1],vals[2],vals[3],vals[4]))
                vetoVetoExpnum.append(vals[0])
                vetoVetoCcd.append(vals[1])
                statvetoVetoA.append(vals[2])
                statvetoVetoB.append(vals[3])
                niteVetoVeto.append(vals[4])
            
        vetoGoodExpnum = []
        vetoGoodCcd = []
        statvetoGoodA = []
        statvetoGoodB = []
        niteVetoGood = []

        with open('vetogood.txt', 'w') as veto:
            for vals in vetoGood:
                veto.write("%10d %2d %10s %7s %7s \n" % (vals[0],vals[1],vals[2],vals[3],vals[4]))
                vetoGoodExpnum.append(vals[0])
                vetoGoodCcd.append(vals[1])
                statvetoGoodA.append(vals[2])
                statvetoGoodB.append(vals[3])
                niteVetoGood.append(vals[4])

        with open('vetogood.txt', 'a') as veto:
            for vals in goodVeto:
                veto.write("%10d %2d %10s %7s %7s \n" % (vals[0],vals[1],vals[2],vals[3],vals[4]))
                vetoGoodExpnum.append(vals[0])
                vetoGoodCcd.append(vals[1])
                statvetoGoodA.append(vals[2])
                statvetoGoodB.append(vals[3])
                niteVetoGood.append(vals[4])
            

        #count how many ccd in each exposure are declared veto veto (both amps)
        singleExps = sorted(set(vetoVetoExpnum))

        
        #loop through singleExps and count how many ccds are veto veto for each exposure.
        for exp in singleExps:
            count = 0
            #myccd = []
            for i,vetoexp in enumerate(vetoVetoExpnum):
                if vetoexp == exp:
                    count += 1
                    #myccd.append(vetoVetoCcd[i])
                    #print "ccd = %d", vetoVetoCcd[i]
                else:
                    continue
            #ccd 45 has bright bad columns that affects the rms calculated by Robert
            #If there is an exposure with VETO VETO and the ccd=45, then this is not to be 
            #declared bad exposure.
            #print myccd
            #theccds = ','.join(myccd)
            if vetoVetoCcd[i] == 45 and count == 1:
                #update the database BIAS_STATUS as GOOD
                #update_bias_query = """update ricardoc.bias_qa_v2 SET BIAS_STATUS = 'GOOD' where expnum = %d """ % (exp)
                insert_bias_good_veto = """ insert into ricardoc.bias_good_veto (EXPNUM, BIAS_STATUS) values (%d, 'GOOD') """ % (exp)
                #flag_change_to_good = True
                #flag_change_to_veto = False
                #print "update bias query", update_bias_query
                print "exp %d has %d ccds , ccd= %d as VETO VETO" % (exp, count, vetoVetoCcd[i])
            #if there is more than 1 bad ccd with VETO VETO, then declare the whole exposure as VETO
            elif count > 1:
                insert_bias_good_veto = """insert into ricardoc.bias_good_veto (EXPNUM, BIAS_STATUS) values (%d, 'VETO')""" % (exp)
                #update_bias_query = """update ricardoc.bias_qa_v2 SET BIAS_STATUS = 'VETO' where expnum = %d""" % (exp)
                #flag_change_to_veto = True
                #flag_change_to_good = False
                #print "update VETO VETO to VETO", update_bias_query
                print "exp %d has %d ccds as VETO VETO" % (exp, count) 
            
            insert2Db(dbh, insert_bias_good_veto,verbose=True)
            
            
            #if flag_change_to_good:
            #    for i in myccd:
            #        update_bias_query = """update ricardoc.bias_qa_v2 SET BIAS_STATUS = 'GOOD' where expnum = %d and ccd = %d """ % (exp,i)
            #        print update_bias_query
            #        insert2Db(dbh, update_bias_query,verbose=True)
            #elif flag_change_to_veto:
            #    for i in myccd:
            #        update_bias_query = """update ricardoc.bias_qa_v2 SET BIAS_STATUS = 'VETO' where expnum = %d and ccd = %d""" % (exp,i)
            #        print update_bias_query
            #        insert2Db(dbh, update_bias_query,verbose=True)
            #"""
         
        
        #vetoGoodQuery = """update ricardoc.bias_qa_v2 SET BIAS_STATUS = 'GOOD' where BIAS_STATUS is NULL or bias_status != 'VETO' and expnum in (select distinct(expnum) from bias_qa_v2 where 
        #            bias_status_a = 'VETO' and bias_status_b = 'GOOD' and ccd != 2)"""
        
        vetoGoodQuery = """select distinct(b.expnum) from ricardoc.bias_qa_v2 b, exposure e where 
                        b.bias_status_a = 'VETO' and b.bias_status_b = 'GOOD' 
                        and b.ccd != 2 and e.expnum = b.expnum and e.nite > 20130801 and e.nite < 201402010 and
                        b.expnum not in (select expnum from ricardoc.bias_good_veto where bias_status is not NULL)""" 

        vetogood = query_to_cur(dbh, vetoGoodQuery, verbose=True)
        for exposure in vetogood:
            insert_veto_good = """insert into ricardoc.bias_good_veto (EXPNUM, BIAS_STATUS) VALUES (%d, 'GOOD')""" % (exposure)
            insert2Db(dbh, insert_veto_good, verbose=True)


        goodVetoQuery = """select distinct(b.expnum) from ricardoc.bias_qa_v2 b, exposure e where 
                        b.bias_status_a = 'GOOD' and b.bias_status_b = 'VETO' 
                        and b.ccd != 2 and e.expnum = b.expnum and e.nite > 20130801 and e.nite < 201402010
                        and b.expnum not in (select expnum from ricardoc.bias_good_veto where bias_status is not NULL)""" 

        goodveto = query_to_cur(dbh, goodVetoQuery)
        for exposure in goodveto:
            insert_good_veto = """insert into ricardoc.bias_good_veto (EXPNUM, BIAS_STATUS) VALUES (%d, 'GOOD')""" % (exposure)
            insert2Db(dbh, insert_good_veto, verbose=True)

        goodGoodQuery = """select distinct(b.expnum) from ricardoc.bias_qa_v2 b, exposure e where 
                        b.bias_status_a = 'GOOD' and b.bias_status_b = 'GOOD' 
                        and b.ccd != 2 and e.expnum = b.expnum and e.nite > 20130801 and e.nite < 201402010 and
                        b.expnum not in (select expnum from ricardoc.bias_good_veto where bias_status is not NULL)""" 

        goodgood = query_to_cur(dbh, goodGoodQuery, verbose=True)
        for exposure in goodgood:
            insert_good_good = """insert into ricardoc.bias_good_veto (EXPNUM, BIAS_STATUS) VALUES (%d, 'GOOD')""" % (exposure)
            insert2Db(dbh, insert_good_good, verbose=True)


        #goodVetoQuery = """update ricardoc.bias_qa_v2 SET BIAS_STATUS = 'GOOD' where BIAS_STATUS is NULL or bias_status != 'VETO' and expnum in (select distinct(expnum) from bias_qa_v2  where 
        #            bias_status_a = 'GOOD' and bias_status_b = 'VETO' and ccd != 2)"""

        #goodVetoQuery = """insert into ricardoc.bias_good_veto (BIAS_STATUS) values ('GOOD') where bias_status is NULL and expnum in (select distinct(expnum) from ricardoc.bias_qa_v2  where 
        #                bias_status_a = 'GOOD' and bias_status_b = 'VETO' and ccd != 2)"""

        ##insert2Db(dbh, vetoGoodQuery,verbose=True)
        #insert2Db(dbh, goodVetoQuery,verbose=True)

        #dbh.commit()

        
        print "total Number of exposures with at least one ccd VETO", len(singleExps)
        print "total Number of exposures: ", len(expnum)

    
    #create a list of veto flats and veto bias
    select_veto_flats = """select distinct(expnum) from flats_stats_qa where flat_status='VETO'"""
    select_veto_bias = """select expnum from bias_good_veto where bias_status='VETO'"""
    
    veto_flats = query_to_cur(dbh, select_veto_flats, verbose)
    veto_bias = query_to_cur(dbh, select_veto_bias, verbose)
    
    with open('veto_exposures.list','w') as f:
        for exp in veto_flats:
            f.write("%d \n" % (exp))
        for exp in veto_bias:
            f.write("%d \n" % (exp))

    #create a list of veto flats and veto bias
    select_good_flats = """select distinct(expnum) from flats_stats_qa where flat_status='GOOD'"""
    select_good_bias = """select expnum from bias_good_veto where bias_status='GOOD'"""
    
    good_flats = query_to_cur(dbh, select_good_flats, verbose)
    good_bias = query_to_cur(dbh, select_good_bias, verbose)
    
    with open('good_exposures.list','w') as f:
        for exp in good_flats:
            f.write("%d \n" % (exp))
        for exp in good_bias:
            f.write("%d \n" % (exp))

            
    
    
    print "fisnihed plotting ...."


    """
    create table ricardoc.bias_qa (
    ID INTEGER NOT NULL,
    EXPNUM NUMBER(10) NOT NULL,
    CCD NUMBER(3),
    BIAS_STATUS_A VARCHAR2(25),
    BIAS_STATUS_B VARCHAR2(25),
    BIAS_STATUS VARCHAR2(25),
    constraint bias_qa_pk PRIMARY KEY (ID)
    );
    """

if  __name__ == '__main__':
    main()
