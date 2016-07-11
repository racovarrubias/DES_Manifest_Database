#!/usr/bin/env python
"""Calculates Gain and rdnoise for DES images

.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>

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
parser.add_argument("-ff", "--flat", required=False, dest="flat", help="Input two raw flats, comma separated")
parser.add_argument("-bb", "--bias", required=False, dest="bias", help="Input two raw bias, comma separated")
parser.add_argument("-n", "--nite", required=True, dest="night", help="Night to select the bias and flats automatically. This will override input flats and biases")
parser.add_argument("-f", "--filter", required=False, dest="filter", help="Filter to select flats")
parser.add_argument("-gr", "--grfile", required=False, dest="grfile", default=None,help="Gain_rdoise file created by this code")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")


args = parser.parse_args()
flats = args.flat
biases = args.bias
night = args.night
filter = args.filter
grfile = args.grfile
verbose = args.verbose


def robust_mean(x):
    y = x.flatten()
    if len(y) < 6:
        return -999.
    if len(np.unique(y))==1:
        meany=y[0]
    else:
        n = len(y)
        y.sort()
        ind_qt1 = round((n+1)/4.)
        ind_qt3 = round((n+1)*3/4.)
        IQR = y[ind_qt3]- y[ind_qt1]
        lowFense = y[ind_qt1] - 1.*IQR
        highFense = y[ind_qt3] + 1.*IQR
        if lowFense == highFense:
            meany=lowFense
        else:
            ok = (y>lowFense)*(y<highFense)
            yy=y[ok]
            meany=yy.mean(dtype='double')
    return meany

def robust_std(x):
    y = x.flatten()
    if len(y) < 6:
        return -999
    if len(np.unique(y))==1:
        meany=y[0]
    else:
        n = len(y)
        y.sort()
        ind_qt1 = round((n+1)/4.)
        ind_qt3 = round((n+1)*3/4.)
        IQR = y[ind_qt3]- y[ind_qt1]
        lowFense = y[ind_qt1] - 1.5*IQR
        highFense = y[ind_qt3] + 1.5*IQR
        if lowFense == highFense:
            meany=lowFense
        else:
            ok = (y>lowFense)*(y<highFense)
            yy=y[ok]
            meany=yy.mean(dtype='double')
            stdy = yy.std(dtype='double')
    return stdy

def get_col_row_range(trimseca,trimsecb):
    """Gets the minimum and maximum columns and rows to use for statistics
    
        :param trimseca:    trim section from Amp A from image header
        :param trimsecb:    trim section for Amp B from image header
        :return: list [mincol_A,maxcol_A,minrow_A,maxrow_A,mincol_B,maxcol_B,minrow_B,maxrow_B]
    """
    
    #One example of regular expression to separate the string in a patter x=[X1:X2,X3,X4]
    #m=re.search('\[([^:]+):([^,]+),([^:]+):([^:]+)\]',x)
    #If we only have numbers instead of strings I can do:
    #m=re.search('(\d+):(\d+),(\d+):(\d+)',x
    
    section_A = re.search('(\d+):(\d+),(\d+):(\d+)',trimseca)
    mincol_A = section_A.group(1)
    maxcol_A = section_A.group(2)
    minrow_A = section_A.group(3)
    maxrow_A = section_A.group(4)
    
    section_B = re.search('(\d+):(\d+),(\d+):(\d+)',trimsecb)
    mincol_B = section_B.group(1)
    maxcol_B = section_B.group(2)
    minrow_B = section_B.group(3)
    maxrow_B = section_B.group(4)
    
    
    limits_A = [mincol_A,maxcol_A,minrow_A,maxrow_A]
    limits_B = [mincol_B,maxcol_B,minrow_B,maxrow_B]

    return limits_A,limits_B

def calc_gain_rdnoise(flat1,flat2,bias1,bias2,row_col_A,row_col_B):
    """calculates gain and rdnoise for each amplifier

        :param flat1:    raw flat used to get gain and rdnoise
        :param flat2:    raw flat used to get gain and rdnoise
        :param bias1:    raw bias used to get gain and rdnoise
        :param bias2:    raw bias used to get gain and rdnoise
        :param row_col_A list with min and max col,row where good data is for Amp A
        :param row_col_B list with min and max col,row where good data is for Amp B
        :return: gain and rdnoise for each amplifier
    """
    
    #calculate gain and rdnoise in slices (boxes) of 200 by 200 size.
    box_size_row = 400
    box_size_col = 200
    flatdiff_A =[]
    biasdiff_A = []
    flatdiff_B = []
    biasdiff_B = []
    mean_flat1_A = []
    mean_flat1_B = []
    mean_flat2_A = []
    mean_flat2_B = []
    mean_bias1_A = []
    mean_bias1_B = []
    mean_bias2_A = []
    mean_bias2_B = []
    std_flatdiff_A = []
    std_flatdiff_B = []
    std_biasdiff_A = []
    std_biasdiff_B = []
    A_box =[]
    B_box=[]
    
    #caclualte the mean of flats and bias 
    
    for col in range(int(row_col_A[0])+200, int(row_col_A[1])-200, box_size_col):
        for row in range(int(row_col_A[2])+200, int(row_col_A[3])-200, box_size_row):
            #calculate the flat and bias difference for each box and determine the mean of all the values
            if (row + box_size_row > int(row_col_A[3])-200) or (col + box_size_col > int(row_col_A[1])-200):
                continue 
            flats_dif  = flat1[row:row + box_size_row,col:col + box_size_col] - flat2[row:row + box_size_row,col:col + box_size_col]
            std_flatdiff = robust_std(flats_dif)
            robust_std_flatdiff = robust_std(flats_dif)
            std_flatdiff_A.append(std_flatdiff)
            bias_dif = bias1[row:row + box_size_row,col:col + box_size_col] - bias2[row:row + box_size_row,col:col + box_size_col]
            std_biasdiff = robust_std(bias_dif)
            std_biasdiff_A.append(std_biasdiff)
            mean_flat1_A.append(robust_mean(flat1[row:row + box_size_row,col:col + box_size_col]))
            mean_flat2_A.append(robust_mean(flat2[row:row + box_size_row,col:col + box_size_col]))
            mean_bias1_A.append(robust_mean(bias1[row:row + box_size_row,col:col + box_size_col]))
            mean_bias2_A.append(robust_mean(bias2[row:row + box_size_row,col:col + box_size_col]))
            A_box.append([row,row + box_size_row,col,col + box_size_col])

    for col in range(int(row_col_B[0])+200, int(row_col_B[1])-200, box_size_col):
        for row in range(int(row_col_B[2])+200, int(row_col_B[3])-200, box_size_row):
            #calculate the flat and bias difference for each box and determine the mean of all the values
            if (row + box_size_row > int(row_col_B[3])-200) or (col + box_size_col > int(row_col_B[1])-200):
                continue 
            flat_dif = flat1[row:row + box_size_row,col:col + box_size_col] - flat2[row:row + box_size_row,col:col + box_size_col]
            std_flatdiff = robust_std(flat_dif)
            std_flatdiff_B.append(std_flatdiff)
            bias_dif = bias1[row:row + box_size_row,col:col + box_size_col] - bias2[row:row + box_size_row,col:col + box_size_col]
            std_biasdiff =robust_std(bias_dif)
            std_biasdiff_B.append(std_biasdiff)
            mean_flat1_B.append(robust_mean(flat1[row:row + box_size_row,col:col + box_size_col]))
            mean_flat2_B.append(robust_mean(flat2[row:row + box_size_row,col:col + box_size_col]))
            mean_bias1_B.append(robust_mean(bias1[row:row + box_size_row,col:col + box_size_col]))
            mean_bias2_B.append(robust_mean(bias2[row:row + box_size_row,col:col + box_size_col]))  
            B_box.append([row,row + box_size_row,col,col + box_size_col])
            
 
    
    #print "ccd = %d,   mean_flat1_A: %10.4f  and   mean_flat1_B %10.4f and mean_flat2_A %10.4f and   mean_flat2_B %10.4f" % (i+1,mean_flat1_A,mean_flat1_B,mean_flat2_A,mean_flat2_A)
    #print "ccd = %d, std_flatdiff_A: %10.4f  and std_biasdiff_A %10.4f  std_flat_diff_B %10.4f and std_biasdiff_B %10.4f" % (i+1,std_flatdiff_A, std_biasdiff_A,std_flatdiff_B,std_biasdiff_B)
    
    #caclulate the mean of the standard deviation for all the sections measured
    std_flatdiff_all_A = scipy.mean(std_flatdiff_A)
    std_flatdiff_all_B = scipy.mean(std_flatdiff_B)
    std_biasdiff_all_A = scipy.mean(std_biasdiff_A)
    std_biasdiff_all_B = scipy.mean(std_biasdiff_B)
    mean_flat1_all_A = scipy.mean(mean_flat1_A)
    mean_flat2_all_A = scipy.mean(mean_flat2_A)
    mean_flat1_all_B = scipy.mean(mean_flat1_B)
    mean_flat2_all_B = scipy.mean(mean_flat2_B)
    mean_bias1_all_A = scipy.mean(mean_bias1_A)
    mean_bias2_all_A = scipy.mean(mean_bias2_A)
    mean_bias1_all_B = scipy.mean(mean_bias1_B)
    mean_bias2_all_B = scipy.mean(mean_bias2_B)
    
    #print np.isinf(mean_flat1_A)
    #print np.isinf(mean_flat1_B)
    #print np.isinf(mean_flat2_A)
    #print np.isinf(mean_flat2_A)
    #print np.isnan(mean_flat1_A)
    #print np.isnan(mean_flat1_B)
    #print np.isnan(mean_flat2_A)
    #print np.isnan(mean_flat2_B)
    #calculate the mean value for all the sections of flats and biases:
    
    gain_A = ((mean_flat1_all_A + mean_flat2_all_A) - (mean_bias1_all_A + mean_bias2_all_A)) / (std_flatdiff_all_A**2 - std_biasdiff_all_A**2)
    gain_B = ((mean_flat1_all_B + mean_flat2_all_B) - (mean_bias1_all_B + mean_bias2_all_B)) / (std_flatdiff_all_B**2 - std_biasdiff_all_B**2)
    
    rdnoise_A = (gain_A * std_biasdiff_all_A) / math.sqrt(2)
    rdnoise_B = (gain_B * std_biasdiff_all_B) / math.sqrt(2)
    
    #print "ccd= %d, gain_A: %5.4f, gain_B: %5.4f " %(i+1,gain_A,gain_B)
    #print "ccd= %d, rdnoise_A: %5.4f, rdnoise_B: %5.4f" % (i+1, rdnoise_A, rdnoise_B)
    
    return [gain_A,gain_B,rdnoise_A,rdnoise_B]


def main():

    #reading gain rdnoise file to calculate the median of each ccd flats and write results in a file
    if grfile is not None:
        
        splitingFile = grfile.split('_')
        theFilter = splitingFile[4].split('.')[0]
        print theFilter
        
        with open(grfile,'rb') as values_data:
            for i in range(4):
                next(values_data)
            reader = csv.DictReader(values_data, delimiter=',')
            allData = [r for r in reader]
    
        print len(allData)
        if len(allData) < 50:
            print "skipping night %s", night
            sys.exit(0)
        mean_flat = []
        #calculate robust mean for each flat
        for i in range(len(allData)):
            mean_flat.append((float(allData[i]['MEAN FLAT1 A']) + float(allData[i]['MEAN FLAT1 B']))/2.0)
        
        meanFlatValuesName = 'mean_flats_' + theFilter + '.txt'
        with open(meanFlatValuesName,'a') as outfile:
            outfile.write("%.1f %.1f %s %s \n" % (np.median(mean_flat), np.std(mean_flat), theFilter, night)) 
        
        
        #print np.mean(mean_flat)
        #print np.std(mean_flat)

        sys.exit(0)

    #print "band:" ,filter
    #directory where to look for raw data
    if inputData = 'src':
        #src data will be used to calcualte gain and rdnoise
        root_dir = '/archive_data/Archive/DTS/src/'
    elif inputData = 'raw':
        #crostalked, overscan subtracted flats and bias will be used to calculate gain and rdnoise
        root_dir = '/archive_data/Archive/OPS/'
        
    if flats is None:
        print "No flats given"
        night_flag = True
    else:
        (flat1,flat2) = flats.split(',')
        print "flat1 : %s and flat2 : %s \n" % (flat1,flat2)
        night_flag = False
        
    if biases is None:
        print "No biases given"
        night_flag = True
    else:
        (bias1,bias2) = biases.split(',')
        print "bias1 : %s and bias2 : %s \n" % (bias1,bias2)
        night_flag = False

    if night_flag == True:
        print "Flats and bias will be automatically selected from the input night. Command-line arguments given for flats and bias will not be used."
    else:
        night_flag = False
    

    
    try:
        desdmfile = os.environ["DES_SERVICES"]
    except KeyError:
        desdmfile = None
    dbh = coreutils.desdbi.DesDbi(desdmfile,"db-desoper")
    cur_flats = dbh.cursor()
    cur_biases = dbh.cursor()


    queryitems = ["id","band","exposurename","exptime","nite","mjd_obs"]
    querylist = ",".join(queryitems)


    #DEfine a fixed band... I will only use r band data
    if filter is None:
        band = 'r'
    else:
        band = filter
    
    if band =='r' or band == 'z' or band == 'Y':
        band_exptime = 10.0
    elif band == 'i':
        band_exptime = 22.0
    elif band == 'u' or band == 'g':
        band_exptime = 30.0
    else:
        print "filter must be r, i, z or Y"
        sys.exit(1) 

    #print band
    
    if night_flag:
        query_flats = """select %s from exposure where 
                        obstype = 'dome flat' and 
                        band = '%s' and 
                        project='DTS' 
                        and nite = '%s' order by exposurename""" % (querylist, band, night)
        
        query_biases = """ select %s from exposure where 
                        obstype = 'zero' and
                        project = 'DTS' and
                        object IS NULL and
                        band like '%%Empty%%' and
                        nite = '%s' order by exposurename""" % (querylist, night)
    #else:
    #    query_flats = """ select %s from exposure
    #                    where  obstype = 'dome flat' and
    #                    exposurename ='%s' and
    #                    project = 'DTS' and
    #                    exposurename ='%s' """ % (querylist, f1, f2)
    #
    #    query_biases = """ select %s from exposure where 
    #                    obstype = 'zero' and
    #                    project = 'DTS' and  
    #                    band like '%%Empty%%' and
    #                    exposurename = '%s' and
    #                    exposurename = '%s' """ % (querylist, b1, b2)

        if verbose:
            print "flats query:\n", query_flats
            print "\n bias query:\n", query_biases
            print "\n executing queries"
        
        #execute the query
        cur_flats.arraysize = 1000 # get 1000 at a time when fetching
        cur_flats.execute(query_flats)
    
        cur_biases.arraysize = 1000
        cur_biases.execute(query_biases)



    if night_flag:
        #grab a couple of flats
        #select two flats. For safety, don;t select first flat from lis.
        flat1 = None
        flat2 = None
        bias1 = None
        bias2 = None
        
        #while (flat1 is None) or (flat2 is None):
        for flat in cur_flats:
            #print flat
            exptime = flat[3]
            if band == 'r' or band == 'z' or band == 'Y':
                if exptime == band_exptime and flat1 is None:
                    flat1 = flat[2]
                elif (exptime == band_exptime) and (flat2 is None):
                    flat2 = flat[2]
                else:
                    continue
            if band == 'i':
                if exptime == band_exptime and flat1 is None:
                    flat1 = flat[2]
                elif (exptime == band_exptime) and (flat2 is None):
                    flat2 = flat[2]
                else:
                    continue
            if band == 'u' or band == 'g':
                if exptime == band_exptime and flat1 is None:
                    flat1 = flat[2]
                elif (exptime == band_exptime) and (flat2 is None):
                    flat2 = flat[2]
                else:
                    continue

    
    
        #select two bias. For safety, don;t select first bias from lis.
        bcount=0
        for bias in cur_biases:
            mjd_obs = bias[5]
            if bcount == 0:
                elast = 0.0
                lastmjd = mjd_obs
            else:
                elast = 3600.*24.*(mjd_obs - lastmjd)
                lastmjd = mjd_obs
            #print lastmjd, elast
            if bcount > 1 and elast < 30. and bias1 is None:
                #if bcount == 2:
                bias1 = bias[2]
                mycount_b1 = bcount
                #print "mycount", mycount_b1
            if bcount > 1 and elast < 30. and bias2 is None:
                if bcount > mycount_b1:
                    #print bcount
                    bias2 = bias[2]
                
            bcount += 1
            if (bias1 is not None) and (bias2 is not None):
                break
    
    print "flats selected: %s and %s" % (flat1,flat2)
    print "bias selected: %s and %s" % (bias1, bias2)
    
    #I am assuming all images are src, from DTS.
    #copy the images to a tmp directory /tmp
    #and read them
    flat1_loc = root_src_dir + night + '/src/' + flat1 + '.fz'
    flat2_loc = root_src_dir + night + '/src/' + flat2 + '.fz'
    bias1_loc = root_src_dir + night + '/src/' + bias1 + '.fz'
    bias2_loc = root_src_dir + night + '/src/' + bias2 + '.fz'

    os.system("cp %s /tmp/" % flat1_loc)
    os.system("cp %s /tmp/" % flat2_loc)
    os.system("cp %s /tmp/" % bias1_loc)
    os.system("cp %s /tmp/" % bias2_loc)
    
    tmp_flat1 = '/tmp/' + flat1 + '.fz'
    tmp_flat2 = '/tmp/' + flat2 + '.fz'
    tmp_bias1 = '/tmp/' + bias1 + '.fz'
    tmp_bias2 = '/tmp/' + bias2 + '.fz'

    hdu_flat1 = pf.open(tmp_flat1)
    hdu_flat2 = pf.open(tmp_flat2)
    hdu_bias1 = pf.open(tmp_bias1)
    hdu_bias2 = pf.open(tmp_bias2)
    
    trimsec_A = hdu_flat1[1].header['TRIMSECA'] 
    trimsec_B = hdu_flat1[1].header['TRIMSECB']
        
    #determine Size of each side of ccd and get numbers of columns and rows.
    (col_row_A,col_row_B) = get_col_row_range(trimsec_A,trimsec_B)

    
    data_flat1 = {}
    data_flat2 = {}
    data_bias1 = {}
    data_bias2 = {}
    gain_A = {}
    gain_B ={}
    rdnoise_A ={}
    rdnoise_B ={}
    mean_flat1_A = {}
    mean_flat1_B = {}
    mean_flat2_A = {}
    mean_flat2_B = {}
    mean_bias1_A = {}
    mean_bias1_B = {}
    mean_bias2_A = {}
    mean_bias2_B = {}
    night_data = {}

    #Also calculates basic statistics from input image (mean and std for a large region for each amp
    
    f = open("gain_rdnoise_%s_%s.csv"% (night,band), "w")
    writer = csv.writer(f, delimiter=',',  quotechar='"', quoting=csv.QUOTE_MINIMAL)
    flats_used = ['#Flats Used For Calculation: ', flat1_loc, flat2_loc]
    bias_used = ['#Bias Used For Calculation: ', bias1_loc, bias2_loc]
    sectionA = ['#Section to calculate mean in Amp A: [1300:1600,1800:2800]']
    sectionB = ['#Section to calculate mean in Amp B: [400:700,1800:2800]']
    writer.writerow(flats_used)
    writer.writerow(bias_used)
    writer.writerow(sectionA)
    writer.writerow(sectionB)
    writer.writerow(['CCD', 'GAIN A', 'GAIN B', 'RDNOISE A', 'RDNOISE B', 'MEAN FLAT1 A', 'MEAN FLAT1 B', 'MEAN FLAT2 A', 'MEAN FLAT2 B', 'MEAN BIAS1 A', 'MEAN BIAS1 B', 'MEAN BIAS2 A', 'MEAN BIAS2 B','NIGHT']) #csv Header

    if verbose:
        print "trimsec A: %s \ntrimsec B: %s" % (trimsec_A, trimsec_B)
        text_sectiona = '[1300:1600,1800:2800]'
        text_sectionb = '[400:700,1800:2800]'
        print "Calculating mean of flats and bias for Amp A and B using section A %s and section B %s\n" % (text_sectiona,text_sectionb)
    
    for i in range(62):
        
        if (i+1) == 61: # or (i+1) == 02:
            continue
        data_flat1[i+1] = hdu_flat1[i+1].data
        data_flat2[i+1] = hdu_flat2[i+1].data
        data_bias1[i+1] = hdu_bias1[i+1].data
        data_bias2[i+1] = hdu_bias2[i+1].data
        
        
        if i == 0:
            if verbose:
                print "CCD, GAIN A, GAIN B, RDNOISE A, RDNOISE B, MEAN FLAT1 A, MEAN FLAT1 B, MEAN FLAT2 A, MEAN FLAT2 B, MEAN BIAS1 A, MEAN BIAS1 B, MEAN BIAS2 A, MEAN BIAS2 B, NIGHT"
        
        #Calculate Gain and Rdnoise
        gain_rdnoise = calc_gain_rdnoise(data_flat1[i+1],data_flat2[i+1], data_bias1[i+1], data_bias2[i+1],col_row_A,col_row_B)

        gain_A[i+1] = gain_rdnoise[0]
        gain_B[i+1] = gain_rdnoise[1]
        rdnoise_A[i+1] = gain_rdnoise[2]
        rdnoise_B[i+1] = gain_rdnoise[3]
        
        #These values are used to print the information in the file. This is not used to caclualte gain and rdnoise
        mean_flat1_A[i+1] = robust_mean(data_flat1[i+1][1300:1600,1800:2800])
        mean_flat2_A[i+1] = robust_mean(data_flat2[i+1][1300:1600,1800:2800])
        mean_flat1_B[i+1] = robust_mean(data_flat1[i+1][400:700,1800:2800])
        mean_flat2_B[i+1] = robust_mean(data_flat2[i+1][400:700,1800:2800])

        mean_bias1_A[i+1] = robust_mean(data_bias1[i+1][1300:1600,1800:2800])
        mean_bias2_A[i+1] = robust_mean(data_bias2[i+1][1300:1600,1800:2800])
        mean_bias1_B[i+1] = robust_mean(data_bias1[i+1][400:700,1800:2800])
        mean_bias2_B[i+1] = robust_mean(data_bias2[i+1][400:700,1800:2800])
        
        night_data[i+1] = night
                
        ccd=i+1
        if verbose:
            print (" {0:5d} {1:.2f} {2:.2f} {3:.2f} {4:.2f} {5:.2f} {6:.2f} {7:.2f} {8:.2f} {9:.2f} {10:.2f} {11:.2f} {12:.2f} {13:10s}").format(ccd,gain_A[i+1], gain_B[i+1],rdnoise_A[i+1],rdnoise_B[i+1],mean_flat1_A[i+1],mean_flat1_B[i+1],mean_flat2_A[i+1],mean_flat2_B[i+1],mean_bias1_A[i+1],mean_bias1_B[i+1],mean_bias2_A[i+1],mean_bias2_B[i+1],night_data[i+1])

        line = [ccd,"{0:0.2f}".format(gain_A[i+1]), "{0:0.2f}".format(gain_B[i+1]),"{0:0.2f}".format(rdnoise_A[i+1]),"{0:0.2f}".format(rdnoise_B[i+1]),"{0:0.2f}".format(mean_flat1_A[i+1]),"{0:0.2f}".format(mean_flat1_B[i+1]),"{0:0.2f}".format(mean_flat2_A[i+1]),"{0:0.2f}".format(mean_flat2_B[i+1]),"{0:0.2f}".format(mean_bias1_A[i+1]),"{0:0.2f}".format(mean_bias1_B[i+1]),"{0:0.2f}".format(mean_bias2_A[i+1]),"{0:0.2f}".format(mean_bias2_B[i+1]),night_data[i+1]]
        #format_line = ["{0:0.2f}".format(i) for i in line]
        writer.writerow(line)

        #save the data in a file
        
    #execute query again and print gain and rdnoise from image header.
    #execute the query
    try:
        desdmfile = os.environ["DES_SERVICES"]
    except KeyError:
        desdmfile = None
    dbh = coreutils.desdbi.DesDbi(desdmfile,"db-desoper")
    cur = dbh.cursor()
 
    query_items = ["imagename","ccd","gaina","gainb","rdnoisea","rdnoiseb","run"]
    new_querylist = ",".join(query_items)

    rootname_flat1 = flat1.split('.')[0]
    query = """select %s from image where 
                imagename like '%s%%' 
                and run in (select run from image where imagename like '%s%%' and nite= '%s') order by ccd""" % (new_querylist, rootname_flat1, rootname_flat1,night)
    if verbose:
        print query

    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)

    f = open("gain_rdnoise_db_image_%s.csv"% (night), "w")
    writer = csv.writer(f, delimiter=',',  quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow (query_items) #csv Header
    
    for item  in cur:
        if verbose:
            print item
        writer.writerow(item)
    
    #now calculate gainA, rdnosieA, gainB, rdnoiseB for all ccds. Ommit calculations for ccd = 61
    #the gain calculation is as follow (need to define a section in the ccd to calculate this)
    #flatdiff = flat1 - flat2
    #biasdiff = bias1 - bias2
    #
    # gain = ( (mean_flat1 _ mean_flat2) - (mean_bias1 + mean_bias2) ) / (stdv_flatdiff^2 - stdv_biasdiff^2)
    #rdnoise = gain * stdv_biasdiff / sqr(2)
    os.system("rm -f %s %s %s %s" % (tmp_flat1, tmp_flat2, tmp_bias1, tmp_bias2))
          
    
    #Determines success
    return 0
    
if  __name__ == '__main__':
    status = main()
    sys.exit(status)
    
