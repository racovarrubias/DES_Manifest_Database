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


verbose=True

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


    """    
     EXPNUM                                NUMBER(10)
     FLAT_STATUS                             VARCHAR2(80)
     SSR_A                                  NUMBER(8,2)
     TOT_RSQ_B                             NUMBER(8,2)
     RED_CHI2_B                             NUMBER(8,2)
     BORDERS_RATIO_B                        NUMBER(8,2)
     STEPS_RATIO_B                             NUMBER(8,2)
     FIT_COEFF_ZERO_B                        NUMBER(8,2)
     FIT_COEFF_ONE_B                        NUMBER(8,2)
     FIT_COEFF_TWO_B                        NUMBER(8,2)
     FIT_COEFF_THREE_B                        NUMBER(8,2)
     FIT_COEFF_FOURTH_B                        NUMBER(8,2)
     MEAN_CCD_B                             NUMBER(8,2)
     SSR_B                                  NUMBER(8,2)
     TOT_RSQ_A                             NUMBER(8,2)
     RED_CHI2_A                             NUMBER(8,2)
     BORDERS_RATIO_A                        NUMBER(8,2)
     STEPS_RATIO_A                             NUMBER(8,2)
     FIT_COEFF_ZERO_A                        NUMBER(8,2)
     FIT_COEFF_ONE_A                        NUMBER(8,2)
     FIT_COEFF_TWO_A                        NUMBER(8,2)
     FIT_COEFF_THREE_A                        NUMBER(8,2)
     FIT_COEFF_FOURTH_A                        NUMBER(8,2)
     MEAN_CCD_A                             NUMBER(8,2)
     FLAT_STATUS_A                             VARCHAR2(80)
     FLAT_STATUS_B                             VARCHAR2(80)
    """

    section = 'db-desoper'
    dbh = connectDB(section)

    #plot TOT_RSQ vs reduced chi sqr
    select_data_a = """select tot_rsq_a, red_chi2_a, borders_ratio_a, ssr_a, 
                        mean_ccd_a, fit_coeff_zero_a, fit_coeff_one_a, fit_coeff_two_a,
                        fit_coeff_three_a, fit_coeff_fourth_a from ricardoc.flats_stats_qa where band='g'"""

    select_data_veto_a = """select tot_rsq_a, red_chi2_a, borders_ratio_a, ssr_a, 
                        mean_ccd_a, fit_coeff_zero_a, fit_coeff_one_a, fit_coeff_two_a,
                        fit_coeff_three_a, fit_coeff_fourth_a from ricardoc.flats_stats_qa where band='g' and flat_status='VETO'"""

    select_data_b = """select tot_rsq_b, red_chi2_b, borders_ratio_b, ssr_b,
                        mean_ccd_b, fit_coeff_zero_b, fit_coeff_one_b, fit_coeff_two_b,
                        fit_coeff_three_b, fit_coeff_fourth_b from ricardoc.flats_stats_qa where band='g'"""

    select_data_veto_b = """select tot_rsq_b, red_chi2_b, borders_ratio_b, ssr_b,
                        mean_ccd_b, fit_coeff_zero_b, fit_coeff_one_b, fit_coeff_two_b,
                        fit_coeff_three_b, fit_coeff_fourth_b from ricardoc.flats_stats_qa where band='g' and flat_status='VETO'"""

    
    data_a = query_to_cur(dbh, select_data_a, verbose=False)
    data_b  = query_to_cur(dbh, select_data_b, verbose=False)
    
    #print rsq_chi2
    
    tot_rsq_a = []
    red_chi2_a = []
    borders_ratio_a = []
    ssr_a = []
    mean_ccd_a = []
    fit_coeff_zero_a = []
    fit_coeff_one_a = []
    fit_coeff_two_a = []
    fit_coeff_three_a = []
    fit_coeff_fourth_a = []
    
    for val in select_data_a:
        tot_rsq_a.append(val[0])
        red_chi2_a.append(val[1])
        borders_ratio_a.append(val[2])
        ssr_a.append(val[3])
        mean_ccd_a.append(val[4])
        fit_coeff_zero_a.append(val[5])
        fit_coeff_one_a.append(val[6])
        fit_coeff_two_a.append(val[7])
        fit_coeff_three_a.append(val[8])
        fit_coeff_fourth_a.append(val[9])

    tot_rsq_b = []
    red_chi2_b = []
    borders_ratio_b = []
    ssr_b = []
    mean_ccd_b = []
    fit_coeff_zero_b = []
    fit_coeff_one_b = []
    fit_coeff_two_b = []
    fit_coeff_three_b = []
    fit_coeff_fourth_b = []
        
    for val in select_data_b:
        tot_rsq_b.append(val[0])
        red_chi2_b.append(val[1])
        borders_ratio_b.append(val[2])
        ssr_b.append(val[3])
        mean_ccd_b.append(val[4])
        fit_coeff_zero_b.append(val[5])
        fit_coeff_one_b.append(val[6])
        fit_coeff_two_b.append(val[7])
        fit_coeff_three_b.append(val[8])
        fit_coeff_fourth_b.append(val[9])
        

    fig=plt.figure()
    for i in range(1,5):
        ax = fig.add_subplot(2,2,i)
        ax.plot()
    
    plt.figure(1)
    plt.clf()
    #plt.xticks(rotation=45)
    plt.grid()
    plt.rc("font", size=9)

    
    plt.plot(rsq, redchi2, 'bo')
    plt.plot(vetorsq, vetoredchi2, 'ro')
    plt.xlabel('Total Squared Root Mean')
    plt.ylabel('Red chi2')
    
    plt.show()
    print "fisniihed plotting ...."
    

if  __name__ == '__main__':
    main()
