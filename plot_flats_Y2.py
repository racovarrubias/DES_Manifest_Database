#!/usr/bin/env python
"""
Plots the mean flat counts/s for exposiure per night.
It reads the flat_stats_qa_y2 table

.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>

  Needs to setup the following before running:
  matplotlib
  
  Uses despyutils. In order to have despyutils in my mac, go to the directory where the ups directory:
  Make sure the ups table is not unseting python or setting packages I don't have installed.
  
  cd /Users/ricardo/DESDM/despyutils/trunk
  setup -r .
  
"""
# imports
import coreutils.desdbi
import argparse
import os
import csv
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import datetime as dt
import matplotlib.dates as mdates
import numpy as np
import re
import sys
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
import math
import matplotlib.patches
from matplotlib.patches     import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle
import matplotlib.colorbar as cbar


#Fontsize for plots
rcParams['axes.labelsize'] = 7
rcParams['xtick.labelsize'] = 7
rcParams['ytick.labelsize'] = 7
rcParams['legend.fontsize'] = 7

#Typeface for plots
#rcParams['font.family'] = 'serif'
#rcParams['font.serif'] = ['Computer Modern Roman']
#rcParams['text.usetex'] = True

#Tick locator
my_locator = MaxNLocator(6)

# Set up axes and plot some awesome science
#ax.yaxis.set_major_locator(my_locator)


parser = argparse.ArgumentParser(description="Plot mean counts of flats in a range of dates with values normalized to 20140820")
parser.add_argument("-d", "--RangeDates", required=True, dest="rangedates", help="Range of dates to plot. Input ranges as date1-date2. If single date is given, first Date is 20140817 by default")
parser.add_argument("-c", "--ccdnum", required=False, dest="ccdnum", help="CCD to process")
parser.add_argument("-b", "--band", required=True, dest="band", help="band to use. A list of band comma separated is accepted")
parser.add_argument("-v", "--verbose", required=False,action="store_true", help="Use --verbose for verbose output")

args = parser.parse_args()
firstDate,lastDate = args.rangedates.split('-')
band = args.band
ccdnum = args.ccdnum
verbose = args.verbose


################################################ Hard Coded constants
### for DECam, a set of list/arrays and dictionaries are hard-coded
### here
#
#
# CCDSECTIONS:
# A dictionary of all of the CCDSECTIONS from a DECam image, produced like:
# This is how I got the CCDSECTIONS for an exposure:
#  gethead -up DETSEC CCDNUM CRPIX1 CRPIX2 CRVAL1 CRVAL2 PIXSCAL1 DE*[0-9].fits
# | awk '{printf \"%s    %-30s %2d    %12.2f %12.2f %14.8f %14.8f\n\",$1,$2,$3,$4,$5,$6,$7}'
#

CCDSECTIONS = {
    1  : [2049,4096,8193,12288],
    2  : [2049,4096,12289,16384],
    3  : [2049,4096,16385,20480],
    4  : [4097,6144,6145,10240],
    5  : [4097,6144,10241,14336],
    6  : [4097,6144,14337,18432],
    7  : [4097,6144,18433,22528],
    8  : [6145,8192,4097,8192],
    9  : [6145,8192,8193,12288],
    10 : [6145,8192,12289,16384],
    11 : [6145,8192,16385,20480],
    12 : [6145,8192,20481,24576],
    13 : [8193,10240,2049,6144],
    14 : [8193,10240,6145,10240],
    15 : [8193,10240,10241,14336],
    16 : [8193,10240,14337,18432],
    17 : [8193,10240,18433,22528],
    18 : [8193,10240,22529,26624],
    19 : [10241,12288,2049,6144],
    20 : [10241,12288,6145,10240],
    21 : [10241,12288,10241,14336],
    22 : [10241,12288,14337,18432],
    23 : [10241,12288,18433,22528],
    24 : [10241,12288,22529,26624],
    25 : [12289,14336,1,4096],
    26 : [12289,14336,4097,8192],
    27 : [12289,14336,8193,12288],
    28 : [12289,14336,12289,16384],
    29 : [12289,14336,16385,20480],
    30 : [12289,14336,20481,24576],
    31 : [12289,14336,24577,28672],
    32 : [14337,16384,1,4096],
    33 : [14337,16384,4097,8192],
    34 : [14337,16384,8193,12288],
    35 : [14337,16384,12289,16384],
    36 : [14337,16384,16385,20480],
    37 : [14337,16384,20481,24576],
    38 : [14337,16384,24577,28672],
    39 : [16385,18432,2049,6144],
    40 : [16385,18432,6145,10240],
    41 : [16385,18432,10240,14335],
    42 : [16385,18432,14336,18431],
    43 : [16385,18432,18432,22527],
    44 : [16385,18432,22528,26623],
    45 : [18433,20480,2049,6144],
    46 : [18433,20480,6145,10240],
    47 : [18433,20480,10240,14335],
    48 : [18433,20480,14336,18431],
    49 : [18433,20480,18432,22527],
    50 : [18433,20480,22528,26623],
    51 : [20481,22528,4097,8192],
    52 : [20481,22528,8193,12288],
    53 : [20481,22528,12289,16384],
    54 : [20481,22528,16385,20480],
    55 : [20481,22528,20481,24576],
    56 : [22529,24576,6145,10240],
    57 : [22529,24576,10240,14335],
    58 : [22529,24576,14336,18431],
    59 : [22529,24576,18432,22527],
    60 : [24577,26624,8193,12288],
    61 : [24577,26624,12289,16384],
    62 : [24577,26624,16385,20480]
    }

CCDSECTION_X0 = (CCDSECTIONS[28][1]+CCDSECTIONS[35][0])/2.0
CCDSECTION_Y0 = (CCDSECTIONS[35][2]+CCDSECTIONS[28][3])/2.0

# Create the TRIM_CCDSECTION
TRIM_CCDSECTIONS = CCDSECTIONS.copy()
borderpix = 104 # 208/2. as 208 is the space between chips in pixels
for _k,_v in TRIM_CCDSECTIONS.items():
    (_x1,_x2,_y1,_y2) = _v
    _x1 = _x1 + borderpix
    _x2 = _x2 - borderpix
    _y1 = _y1 + borderpix
    _y2 = _y2 - borderpix
    TRIM_CCDSECTIONS[_k] = [_x1,_x2,_y1,_y2]

# Create the Corners of DECam Footprint
DECAM_CORNERS_X = []
DECAM_CORNERS_Y = []
# Bottom [L/R]
for _k in [1,4,8,13,19,25,32,39,45,51,56,60]:
    DECAM_CORNERS_X.append(CCDSECTIONS[_k][0])
    DECAM_CORNERS_X.append(CCDSECTIONS[_k][1])
    DECAM_CORNERS_Y.append(CCDSECTIONS[_k][2])
    DECAM_CORNERS_Y.append(CCDSECTIONS[_k][2])
# Top [L/R]
for _k in [62,59,55,50,44,38,31,24,18,12,7,3]:
    DECAM_CORNERS_X.append(CCDSECTIONS[_k][1])
    DECAM_CORNERS_X.append(CCDSECTIONS[_k][0])
    DECAM_CORNERS_Y.append(CCDSECTIONS[_k][3])
    DECAM_CORNERS_Y.append(CCDSECTIONS[_k][3])
# close the points
DECAM_CORNERS_X.append(DECAM_CORNERS_X[0])
DECAM_CORNERS_Y.append(DECAM_CORNERS_Y[0])
# Into numpy arrays
DECAM_CORNERS_X = np.array(DECAM_CORNERS_X)
DECAM_CORNERS_Y = np.array(DECAM_CORNERS_Y)

######  End of hard-coded constants for DECam #######


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

def drawDECamCCDs_Plot(x0,y0,trim=True,rotate=True,label=False, norm=None, exps=None,**kwargs):

    """ Draws DECam CCDs shapes using matplotlib Plot function on the current plot"""
    #ax = plt.gca()
    fig, ax = plt.subplots()
    allcolors = []
    if trim:
        SECTIONS = TRIM_CCDSECTIONS
    else:
        SECTIONS = CCDSECTIONS
    for k,v in SECTIONS.items():
        (x1,x2,y1,y2) = v
        if rotate:
            x1,y1 = rotate_xy(x1,y1,theta=-90,x0=x0,y0=y0)
            x2,y2 = rotate_xy(x2,y2,theta=-90,x0=x0,y0=y0)

            delta_y = (y2-y1)/2
            delta_x = (x2-x1)/19
            x1bin = []
            x2bin = []
            deltay_A = []
            deltay_B = []
            deltax = []
            
            y1bin = []
            y2bin = []
            colors = []
            
            #ampA
            for bin in range(19):
                x1bin.append(x1 + bin*delta_x)
                x2bin.append(x1+delta_x + bin*delta_x)
                deltax.append(delta_x)
                y1bin.append(y1)
                y2bin.append(y1+delta_y)
                deltay.append(delta_y)
                if k == 2 or k == 61:
                    continue
                colors.append(norm[k][bin])
            
            #AmpB
            for bin in range(19,38):
                x1bin.append(x1 + bin*delta_x)
                x2bin.append(x1+delta_x + bin*delta_x)
                deltax.append(delta_x)
                y1bin.append(y1)
                y2bin.append(y1+delta_y)
                deltay.appedn(delta_y)
                if k == 2 or k == 61:
                    continue
                colors.append(norm[k][bin])
            
            
            #make an arrat with same number of bins for y
        else:
            x1,y1 = rotate_xy(x1,y1,theta=0,x0=x0,y0=y0)
            x2,y2 = rotate_xy(x2,y2,theta=0,x0=x0,y0=y0)

        
        # Into numpy arrays
        x = np.array([x1,x2,x2,x1,x1])
        y = np.array([y1,y1,y2,y2,y1])
        
        if k == 2 or k == 61:
            ax.plot(x,y,**kwargs)
            ax.text(0.5*(x2+x1),0.5*(y2+y1),"%s" % k, ha='center',va='center',fontsize=6)
            continue
        
        colors = np.asarray(colors)
        normal = plt.Normalize(colors.min(), colors.max())
        mycolors = plt.cm.cool(normal(colors))

        allcolors.append(colors)

        ax.plot(x,y,**kwargs)
        for xb,yb,w,h,c in zip(x1bin,y1bin,deltax,deltay_A,mycolors):
            rect = plt.Rectangle((xb,yb),w,h,color=c)
            ax.add_patch(rect)

        if label:
            ax.text(0.5*(x2+x1),0.5*(y2+y1),"%s" % k, ha='center',va='center',fontsize=6)
        
   
    
    singleListValues = np.asarray(reduce(lambda x,y: x+y,allcolors))
    norm = plt.Normalize(singleListValues.min(), singleListValues.max())
   
    fig.colorbar(rect,ax=ax)
    #plt.xlim(-15000,15000)
    #plt.ylim(-15000,15000)
    #plt.xlabel('x-pixels')
    #plt.ylabel('y-pixels')
    #plt.title(exp)
    #plt.colorbar(norm,orientation='vertical')
    #cax, _ = cbar.make_axes(ax) 
    #cb2 = cbar.ColorbarBase(cax, cmap=plt.cm.cool,norm=norm) 

    plt.show()


    return norm

def drawDECamCorners_Plot(x0,y0,rotate=True,**kwargs):

    """ Draws DECam Corners using matplotlib Plot function on the current plot"""

    ax = plt.gca()
    if rotate:
        x,y = rotate_xy(DECAM_CORNERS_X,DECAM_CORNERS_Y,theta=-90,x0=x0,y0=y0)
    else:
        x,y = rotate_xy(DECAM_CORNERS_X,DECAM_CORNERS_Y,theta=  0,x0=x0,y0=y0)
    ax.plot(x,y,**kwargs)
    return


def rotate_xy(x,y,theta,x0=0,y0=0,units='degrees'):
    """
    Rotates (x,y) by angle theta and (x0,y0) translation
    """
    if units == 'degrees':
        d2r = math.pi/180. # degrees to radians shorthand
        theta = d2r*theta
    x_new =  (x-x0)*math.cos(theta) - (y-y0)*math.sin(theta)
    y_new =  (x-x0)*math.sin(theta) + (y-y0)*math.cos(theta)
    return x_new,y_new


def referenceNight(dbh,firstDate, ccd, band):
    """Use first date as the zeropoint for all calculations"""
    
    valuesDict = {}
    
    queryReference = """select  fs.filename,fs.ccdnum,fs.band,fs.nite,tp.domeflor,tp.zd,tp.DOMEAZ,tp.AZ,tp.DOMELOW,tp.DOMEHIGH, fs.mean_ccd_a, fs.mean_ccd_b, (fs.mean_ccd_a + fs.mean_ccd_b)/2
                    from exposure_tel_qa_y2 tp, flats_stats_qa_y2 fs 
                    where fs.nite = '%s' 
                    and fs.ccdnum= %s
                    and fs.band= '%s'
                    and fs.flat_status='GOOD'
                    and fs.filename=tp.filename""" % (firstDate, ccd, band)
                
    cur_reference = query_to_cur(dbh, queryReference, verbose)
    
    for values in cur_reference:
        filename = values[0]
        digit = re.match('DECam_(\d+).fits',filename)
        expnum = int(digit.group(1))
        
        valuesDict[expnum] = values

    return valuesDict
    

def collectData(dbh,firstDate,lastDate,ccd,band):
    """
    Select Mean for each ccd/amp (Mean_CCD_A, MEAN_CCD_B) from flats_stats_qa_y2
    """
    
    valuesDict = {}
                    
    queryValues = """select fs.filename,fs.ccdnum,fs.band,fs.nite,tp.domeflor,tp.ZD,tp.DOMEAZ,tp.AZ,tp.DOMELOW,tp.DOMEHIGH, fs.mean_ccd_a, fs.mean_ccd_b, (fs.mean_ccd_a + fs.mean_ccd_b)/2
                    from exposure_tel_qa_y2 tp, flats_stats_qa_y2 fs 
                    where fs.nite >= '%s'
                    and fs.nite <= '%s'
                    and fs.ccdnum= %s
                    and fs.band= '%s'
                    and fs.flat_status='GOOD'
                    and fs.filename=tp.filename""" % (firstDate, lastDate, ccd, band)
                
    cur_values = query_to_cur(dbh, queryValues, verbose)
    
    for values in cur_values:
        filename = values[0]
        digit = re.match('DECam_(\d+).fits',filename)
        expnum = int(digit.group(1))
        
        valuesDict[expnum] = values
    
    return valuesDict

def collectAllCcdsData(dbh,firstDate,lastDate,band):
    """
    Select Mean for each ccd/amp (Mean_CCD_A, MEAN_CCD_B) from flats_stats_qa_y2
    """
    
    #This dictionary contains all the exposures with their respective values for all ccds, and regions.
    valuesExposures = {}
    valsCcds = {}
    medianExposure = {}
    normRegionsExposures = {}
    normRegionsCcds = {}
    
    medianBoxes = []
                    
    queryValues = """select filename,ccdnum,band,nite, MEDIAN_GBOX1_A, MEDIAN_GBOX2_A, MEDIAN_GBOX3_A,MEDIAN_GBOX4_A, MEDIAN_GBOX5_A, MEDIAN_GBOX6_A, MEDIAN_GBOX7_A,
                    MEDIAN_GBOX8_A, MEDIAN_GBOX9_A, MEDIAN_GBOX10_A, MEDIAN_GBOX11_A, MEDIAN_GBOX12_A, MEDIAN_GBOX13_A, MEDIAN_GBOX14_A,MEDIAN_GBOX15_A, MEDIAN_GBOX16_A, 
                    MEDIAN_GBOX17_A, MEDIAN_GBOX18_A, MEDIAN_GBOX19_A, 
                    MEDIAN_GBOX1_B, MEDIAN_GBOX2_B, MEDIAN_GBOX3_B,MEDIAN_GBOX4_B, MEDIAN_GBOX5_B, MEDIAN_GBOX6_B, MEDIAN_GBOX7_B,
                    MEDIAN_GBOX8_B, MEDIAN_GBOX9_B, MEDIAN_GBOX10_B, MEDIAN_GBOX11_B, MEDIAN_GBOX12_B, MEDIAN_GBOX13_B, MEDIAN_GBOX14_B,MEDIAN_GBOX15_B, MEDIAN_GBOX16_B, 
                    MEDIAN_GBOX17_B, MEDIAN_GBOX18_B, MEDIAN_GBOX19_B
                    from flats_stats_qa_y2
                    where nite >= '%s'
                    and nite <= '%s'
                    and band= '%s'
                    and flat_status='GOOD'
                    order by filename, ccdnum""" % (firstDate, lastDate, band)
                
    cur_values = query_to_cur(dbh, queryValues, verbose)
    
    #Calculates the mean for each exposure using all the regions for all the ccds
    for values in cur_values:
        filename = values[0]
        digit = re.match('DECam_(\d+).fits',filename)
        expnum = int(digit.group(1))
        ccdnum = values[1]
        
        valsCcds[ccdnum] = values[2:]
        valuesExposures[expnum] = valsCcds
        medianBoxes.append(np.asarray(values[4:]))
        if ccdnum == 62:
            mexposure = np.median(medianBoxes)
            print expnum,mexposure
            medianExposure[expnum] = mexposure

    #Normalized all regions for each ccd with respect to its exposure            
    for exp,values in valuesExposures.items():
        for ccd, gboxes in values.items():
            #print ccd, np.asarray(gboxes[2:])/medianExposure[exp]
            norm = np.asarray(gboxes[2:])/medianExposure[exp]
            transformToList = norm.tolist()
            normRegionsCcds[ccd] = transformToList 
        normRegionsExposures[exp] = normRegionsCcds
    
    #normRegionsExposures contains for each exposure all ccds normalized to their respective exposure    
    #print normRegionsExposures.keys()
    #print normRegionsExposures.values()
        
        
    return normRegionsExposures


def calculateRatioRegionExposure(allNightsDataValues,valuesKeys):
    """Calculate the mean for the whole reference exposure.
    For every exposure, calculate the mean for all ccds, me_j
    For each region c_i, calculate the ratio r_ij between the region and
    the mean exposure r_ij = c_i/me_j
    allNightsDataValues[ccdnum][expnums][valuesKeys['MEAN_A']]
    For each exposure need to get the mean counts.
    Then loop through all regions for each ccd/exposure.
    """
    
    
    
    
    

def main():

    #Plot all camera
    # These are the predefined centers of DECam
    x0 = CCDSECTION_X0
    y0 = CCDSECTION_Y0
    #drawDECamCCDs_Plot(x0,y0,rotate=True,label=True,color='black',lw=0.5,ls='-')
    #sys.exit()

    #Connect to operations DB
    dbh = connectDB('db-desoper')

    #If a list of comma separated bands is given, split and create a list
    if len(band.split(',')) > 1:
        filter = band.split(',')
        flagBands = True
        if verbose:
            print "bands", filter
    else:
        filter = band
        flagsBands = False
        if verbose:
            print "only One band provided", band

    allccds = range(1,63)
    #remove ccd 2 and 61 from list of ccds
    allccds.remove(2)
    allccds.remove(61)
    #print allccds
    #Get mean for all ccd in the daterange specified
    #select fs.filename,fs.ccdnum,fs.band,fs.nite,tp.domeflor,tp.ZD,tp.DOMEAZ,tp.AZ,tp.DOMELOW,tp.DOMEHIGH, fs.mean_ccd_a, fs.mean_ccd_b, (fs.mean_ccd_a + fs.mean_ccd_b)/2
    #keys = ['FILENAME','EXPNUM','CCDNUM','BAND','NITE','DOMEFLOOR','ZD','DOMEAZ','AZ','DOMELOW','DOMEHIGH', 'MEAN_A', 'MEAN_B', 'MEAN_AB']
    #Dictionary with names for each position in list of values
    valuesKeys = { 'FILENAME' : 0,
                   'CCDNUM' : 1,
                   'BAND' : 2,
                   'NITE' : 3,
                   'DOMEFLOOR' : 4,
                   'ZD' : 5,
                   'DOMEAZ' : 6,
                   'AZ' : 7,
                   'DOMELOW' : 8,
                   'DOMEHIGH' : 9,
                   'MEAN_A' : 10,
                   'MEAN_B' : 11,
                   'MEAN_AB' : 12 }
    
    
    #dictionary with all values for each exposure, ccd, band
    allNormExposuresCcdsRegions = collectAllCcdsData(dbh,firstDate,lastDate, filter)

    #loop through the exposures and plot each of them
    for exp, normsccds in allNormExposuresCcdsRegions.items():
        print "# Fig 1: DECam corners and CCDs using Polygons "
        #plt.figure(1,figsize=(12,12))
        #ax = plt.subplot(111)
        norm = drawDECamCCDs_Plot(x0,y0,rotate=True,label=True, norm=normsccds, exps=exp, lw=0.5,ls='-')
        plt.xlim(-15000,15000)
        plt.ylim(-15000,15000)
        plt.xlabel('x-pixels')
        plt.ylabel('y-pixels')
        plt.title(exp)
        plt.colorbar(norm,orientation='vertical')
        cax, _ = cbar.make_axes(ax) 
        #cb2 = cbar.ColorbarBase(cax, cmap=plt.cm.cool,norm=norm) 

        plt.show()
        sys.exit()       
        

    #drawDECamCorners_Plot(x0,y0,rotate=True,color='blue',lw=0.5,ls='-')



    referenceDataValues = {}
    refValuesDict = {}
    allNightsDataValues = {}
    dataValuesDict = {}
    #create a dictionary for each
    
    for myband in filter:

        for ccd in allccds:
            #print ccd
            refValuesDict = referenceNight(dbh, firstDate, ccd, myband)
            dataValuesDict = collectData(dbh,firstDate,lastDate, ccd, myband)
            
            referenceDataValues[ccd] = refValuesDict
            allNightsDataValues[ccd] = dataValuesDict
    

        
        #calculate for all ccds the values ratio of each region in amp ccd with respect to mean of exposure
        #It also makes a plot for all data.
        #calculateRatioRegionExposure(allNightsDataValues,valuesKeys)
        
        
        #In order to get extract values from the dictionary, follow this:
        #referenceDataValues[EXPNUM][CCDNUM][Position_in_list_of_data_extracted_from_query]
        #EXPNUM is the key value of the dictionary
        #print referenceDataValues[364478][13][valuesKeys['MEAN_AB']
        
        #keys are all ccds 
        #values is a dictionary where each key is the expnu
        allrefexps = referenceDataValues.values()[0].keys()
        firstExp = min(allrefexps)
        refnight = firstDate
        
        #print referenceDataValues.values()[59].keys()
        #print referenceDataValues.values()[59][349109]
        #print referenceDataValues.keys()
        #print myband
        #rint allrefexps

        #print firstExp
        allexpnums = []
        normcounts = []
        domefloor = []
        tempDomelow = []
        tempDomehi = []
        tempDomefloor = []
        nites = []
        meanA = []
        meanB = []
        
        #ccdnum is the ccd to process
        for expnums in allNightsDataValues.values()[0]:
            allexpnums.append(expnums)
            #print ccdnum,expnums,valuesKeys['MEAN_AB']
            counts = float(allNightsDataValues[int(ccdnum)][expnums][valuesKeys['MEAN_AB']])
            domefloor.append(float(allNightsDataValues[int(ccdnum)][expnums][valuesKeys['DOMEAZ']]))
            tempDomelow.append(float(allNightsDataValues[int(ccdnum)][expnums][valuesKeys['DOMELOW']]))
            tempDomehi.append(float(allNightsDataValues[int(ccdnum)][expnums][valuesKeys['DOMEHIGH']]))
            tempDomefloor.append(float(allNightsDataValues[int(ccdnum)][expnums][valuesKeys['DOMEFLOOR']]))
            nites.append(allNightsDataValues[int(ccdnum)][expnums][valuesKeys['NITE']])
            meanA.append(float(allNightsDataValues[int(ccdnum)][expnums][valuesKeys['MEAN_A']]))
            meanB.append(float(allNightsDataValues[int(ccdnum)][expnums][valuesKeys['MEAN_B']]))
            #print counts
            refcounts = float(referenceDataValues[int(ccdnum)][firstExp][valuesKeys['MEAN_AB']])
            #print refcounts
            normcounts.append(counts/refcounts)
            #print '{0:.2f} {0:.2f} {0:.2f} {0:.2f} {0:10}'.format(expnums,counts,refcounts,counts/refcounts,allNightsDataValues[int(ccdnum)][expnums][valuesKeys['NITE']])
            #if float(counts/refcounts) == 1.0:
                #print "%.2f %.2f %.2f %.4f %s XXXXXXXXX" % (expnums,counts,refcounts,float(counts)/float(refcounts),allNightsDataValues[int(ccdnum)][expnums][valuesKeys['NITE']])
            #else:
            #    print "%.2f %.2f %.2f %.4f %s" % (expnums,counts,refcounts,float(counts)/float(refcounts),allNightsDataValues[int(ccdnum)][expnums][valuesKeys['NITE']])
                
        #Divide each region by the mean of the whole exposure
        
        
        #print allexpnums,len(allexpnums)
        #print normcounts,len(normcounts)
       
        #print firstExp
        #sys.exit()
        #print xdata
        normexps = []
        for val in allexpnums:
            normexps.append(val-firstExp)
    
                        
        expsBigEight = []
        expsLowsix = []
        expsBetween = []
        normBigEight = []
        normLowsix = []
        normBetween = []
        domeBigEight = []
        domeLowsix = []
        domeBetween = []
        azBigEight = []
        azLowsix = []
        azBetween = []
        zdBigEight = []
        zdLowsix = []
        zdBetween = []
        

        #plot exposure versus normcounts as a function of domefloor position
        color = [item*(item/max(domefloor)) for item in domefloor]
        #print 'color',color
        fig, ax = plt.subplots()
        im = ax.scatter(normexps, normcounts, c=color, s=50)
        # Add a colorbar
        fig.colorbar(im, ax=ax)
        #horizontal line
        plt.axhline(y=1.,color='k',ls='dashed')
        #horizontal line
        plt.axhline(y=1.1,color='k',ls='dashed')

        plt.title('Dome Azimuth Position, band %s, ccdnum %s' % (myband, ccdnum))
        mylabel = 'Normalized counts to Exposure %s ' % firstExp
        plt.ylabel(mylabel)
        myxlabel = 'Exposure Number - %s ' % firstExp    
        plt.xlabel(myxlabel)    
        
        # set the color limits - not necessary here, but good to know how.
        #im.set_clim(min(domefloor), max(domefloor))
        domePlotFilename = 'domeAZExpnumCounts_' + myband+ '_'+ ccdnum +'.pdf'
        #plt.show()
        plt.savefig(domePlotFilename)
        plt.close()


        #Plot counts versus dome temperature
        color = [item*(item/max(tempDomehi)) for item in tempDomehi]
        #print 'color',color
        fig, ax = plt.subplots()
        im = ax.scatter(normexps, normcounts, c=color, s=50)
        # Add a colorbar
        fig.colorbar(im, ax=ax)
        #horizontal line
        plt.axhline(y=1.,color='k',ls='dashed')
        #horizontal line
        plt.axhline(y=1.1,color='k',ls='dashed')

        plt.title('Dome Temperature High, band %s, ccdnum %s' % (myband, ccdnum))
        mylabel = 'Normalized counts to Exposure %s ' % firstExp
        plt.ylabel(mylabel)
        myxlabel = 'Exposure Number - %s ' % firstExp    
        plt.xlabel(myxlabel)    
        
        # set the color limits - not necessary here, but good to know how.
        #im.set_clim(min(domefloor), max(domefloor))
        domePlotFilename = 'domeTempHighExpnumCounts_' + myband+'_'+ccdnum +'.pdf'
        #plt.show()
        plt.savefig(domePlotFilename)
        plt.close()

        if myband == 'p':
            #plot one night of flats counts with respect to temperature and see if there is a change.
            #select all counts for nights 20141025, 20141024, 20141023, i band flats
            indx1029 = []
            indx1020 = []
            indx1015 = []
            #print nites
            for i,val in enumerate(nites):
                if val == '20141029':
                    #print i,val
                    indx1029.append(i)
                elif val == '20141020':
                    indx1020.append(i)
                elif val == '20141015':
                    indx1015.append(i)
            
            templow1029 = []
            normcounts1029 = []
            normexps1029 = []
            templow1015 = []
            normcounts1015 = []
            normexps1015 = []
            templow1020 = []
            normcounts1020 = []
            normexps1020 = []
    
            print indx1029
            for val in indx1029:
                templow1029.append(tempDomelow[val])
                normcounts1029.append(normcounts[val])
                normexps1029.append(allexpnums[val])
            for val in indx1015:
                templow1015.append(tempDomelow[val])
                normcounts1015.append(normcounts[val])
                normexps1015.append(normexps[val])
            for val in indx1020:
                templow1020.append(tempDomelow[val])
                normcounts1020.append(normcounts[val])
                normexps1020.append(normexps[val])
    
            
            print "temps 1029", templow1029,normcounts1029,normexps1029
            color1029 = [item*(item/max(templow1029)) for item in templow1029]
            fig, ax = plt.subplots()
            im = ax.scatter(normexps1029, normcounts1029, c=color1029, s=50)
            #plt.xticks([],fontsize=8)
            # Add a colorbar
            cbar = fig.colorbar(im, ax=ax, format="%.2f", spacing='proportional')
            cbar.ax.tick_params(labelsize=10) 
            cbar.set_label('Dome Low Temperature',size=10)
            #calculate the mean and stdev of the counts for the range selected
            ncount = np.asarray(normcounts1029)
            normMean = np.mean(ncount)
            stddev = np.std(ncount)
            #horizontal line
            plt.axhline(y=normMean,color='k',ls='dashed')
            #horizontal line
            plt.axhline(y=normMean+stddev,color='r',ls='dashed')
            plt.axhline(y=normMean-stddev,color='r',ls='dashed')
    
            plt.title('Dome Temperature Low, band %s, ccdnum %s 20141029' % (myband, ccdnum))
            mylabel = 'Normalized counts to Exposure %s ' % firstExp
            plt.ylabel(mylabel)
            myxlabel = 'Exposure Number'  
            plt.xlabel(myxlabel)    
            
            # set the color limits - not necessary here, but good to know how.
            #im.set_clim(min(domefloor), max(domefloor))
            domePlotFilename = 'domeTempLowExpnumCounts_20141029_' + myband+'_'+ccdnum +'.pdf'
            #plt.show()
            plt.savefig(domePlotFilename)
            plt.close()
            
            
            
            #print templow1015,normcounts1015,normexps1015
            color1015 = [item*(item/max(templow1015)) for item in templow1015]
            fig, ax = plt.subplots()
            im = ax.scatter(normexps1015, normcounts1015, c=color1015, s=50)
            # Add a colorbar
            cbar = fig.colorbar(im, ax=ax, format="%.2f",spacing='proportional')
            cbar.ax.tick_params(labelsize=10)
            cbar.set_label('Dome Low Temperature',size=10)
            ncount = np.asarray(normcounts1015)
            normMean = np.mean(ncount)
            stddev = np.std(ncount)
            #horizontal line
            plt.axhline(y=normMean,color='k',ls='dashed')
            #horizontal line
            plt.axhline(y=normMean+stddev,color='r',ls='dashed')
            plt.axhline(y=normMean-stddev,color='r',ls='dashed')
    
            
            #horizontal line
            #plt.axhline(y=1.,color='k',ls='dashed')
            #horizontal line
            #plt.axhline(y=1.1,color='k',ls='dashed')
    
            plt.title('Dome Temperature Low, band %s, ccdnum %s 20141015' % (myband, ccdnum))
            mylabel = 'Normalized counts to Exposure %s ' % firstExp
            plt.ylabel(mylabel)
            myxlabel = 'Exposure Number'    
            plt.xlabel(myxlabel)    
            
            # set the color limits - not necessary here, but good to know how.
            #im.set_clim(min(domefloor), max(domefloor))
            domePlotFilename = 'domeTempLowExpnumCounts_20141015_' + myband+'_'+ccdnum +'.pdf'
            #plt.show()
            plt.savefig(domePlotFilename)
            plt.close()
    
            #print templow1020,normcounts1020,normexps1020
            color1020 = [item*(item/max(templow1020)) for item in templow1020]
            fig, ax = plt.subplots()
            im = ax.scatter(normexps1020, normcounts1020, c=color1020, s=50)
            # Add a colorbar
            cbar = fig.colorbar(im, ax=ax, format="%.2f",spacing='proportional')
            cbar.ax.tick_params(labelsize=10)
            cbar.set_label('Dome Low Temperature',size=10)
            #horizontal line
            plt.axhline(y=1.,color='k',ls='dashed')
            #horizontal line
            plt.axhline(y=1.1,color='k',ls='dashed')
    
            plt.title('Dome Temperature Low, band %s, ccdnum %s 20141020' % (myband, ccdnum))
            mylabel = 'Normalized counts to Exposure %s ' % firstExp
            plt.ylabel(mylabel)
            myxlabel = 'Exposure Number - %s ' % firstExp    
            plt.xlabel(myxlabel)    
            
            # set the color limits - not necessary here, but good to know how.
            #im.set_clim(min(domefloor), max(domefloor))
            domePlotFilename = 'domeTempLowExpnumCounts_20141020_' + myband+'_'+ccdnum +'.pdf'
            #plt.show()
            plt.savefig(domePlotFilename)
            plt.close()
            
        

        #Plot counts versus dome Low temperature
        color = [item*(item/max(tempDomelow)) for item in tempDomelow]
        #print 'color',color
        fig, ax = plt.subplots()
        im = ax.scatter(normexps, normcounts, c=color, s=50)
        # Add a colorbar
        fig.colorbar(im, ax=ax)
        #horizontal line
        plt.axhline(y=1.,color='k',ls='dashed')
        #horizontal line
        plt.axhline(y=1.1,color='k',ls='dashed')

        plt.title('Dome Temperature Low, band %s, ccdnum %s' % (myband, ccdnum))
        mylabel = 'Normalized counts to Exposure %s ' % firstExp
        plt.ylabel(mylabel)
        myxlabel = 'Exposure Number - %s ' % firstExp    
        plt.xlabel(myxlabel)    
        
        # set the color limits - not necessary here, but good to know how.
        #im.set_clim(min(domefloor), max(domefloor))
        domePlotFilename = 'domeTempLowExpnumCounts_' + myband+'_'+ccdnum +'.pdf'
        #plt.show()
        plt.savefig(domePlotFilename)
        plt.close()

        #Plot counts versus dome temperature
        color = [item*(item/max(tempDomefloor)) for item in tempDomefloor]
        #print 'color',color
        fig, ax = plt.subplots()
        im = ax.scatter(normexps, normcounts, c=color, s=50)
        # Add a colorbar
        fig.colorbar(im, ax=ax)
        #horizontal line
        plt.axhline(y=1.,color='k',ls='dashed')
        #horizontal line
        plt.axhline(y=1.1,color='k',ls='dashed')

        plt.title('Dome Temperature Floor, band %s' % myband)
        mylabel = 'Normalized counts to Exposure %s ' % firstExp
        plt.ylabel(mylabel)
        myxlabel = 'Exposure Number - %s ' % firstExp    
        plt.xlabel(myxlabel)    
        
        # set the color limits - not necessary here, but good to know how.
        #im.set_clim(min(domefloor), max(domefloor))
        domePlotFilename = 'domeTempFloorExpnumCounts_' + myband+'_'+ccdnum +'.pdf'
        #plt.show()
        plt.savefig(domePlotFilename)
        plt.close()


        #plt.plot(azBigEight,zdBigEight, 'bo')
        #plt.plot(azLowsix,zdLowsix, 'go')
        #plt.plot(azBetween,zdBetween, 'ro')
        #plt.ylim([49.5,51.5])
        #plt.xlim([0,1])
        #plt.savefig('AzZd.png')
        #plt.close()
        
        mindata = min(normcounts)
        maxdata = max(normcounts)
        #plt.plot(domefloor,normcounts, 'ro')
        plt.ylim([mindata-0.1,maxdata+0.1])
        plt.xlim([min(domefloor)-20,max(domefloor)+20])
        plt.plot(domeLowsix, normLowsix, 'go', markersize=5)
        plt.plot(domeBetween,normBetween,'ro', markersize=5)
        plt.plot(domeBigEight, normBigEight, 'bo', markersize=5)
        #horizontal line
        plt.axhline(y=1.,color='k',ls='dashed')
        plt.axhline(y=1.1,color='k',ls='dashed')
        plt.title('Dome Azimuth Position vs Normalized Counts, band %s, ccdnum %s' % (myband, ccdnum))
        mylabel = 'Normalized counts to Exposure %s ' % firstExp
        plt.ylabel(mylabel)
        myxlabel = 'Dome Azimuth'    
        plt.xlabel(myxlabel)    

        #plt.show()
        filename = 'domeAZVsCounts_'+myband+'_'+ccdnum +'.pdf'
        plt.savefig(filename)
        plt.close()

        

        if myband ==  'g':
            normcounts_g = normcounts
            normexps_g = normexps
        elif myband == 'r':
            normcounts_r = normcounts
            normexps_r = normexps
        elif myband == 'i':
            normexps_i = normexps
            normcounts_i = normcounts
        elif myband == 'z':
            normexps_z = normexps
            normcounts_z = normcounts
        elif myband == 'Y':
            normexps_y = normexps
            normcounts_y = normcounts
    
        #print myband, normexps_g,normcounts_g
    

    for myband in filter:
        #Plot all filters from the range of dates given
        mindata = min(normcounts)
        maxdata = max(normcounts)
    
        #f = plt.figure()
        f, (ax1, ax2, ax3,ax4,ax5) = plt.subplots(5, sharex=True, sharey=False)
        #f.subplots_adjust(hspace=0.001)
        #colorexps = [str(item/max(normexps)) for item in normexps]
        ax1.set_ylim([0.9,1.3])
        #if myband == 'g':
        ax1.set_title('Normalized counts to Exp %s, 20140817, ccdnum %s' % (firstExp, ccdnum))
        ax1.plot(normexps_g,normcounts_g,'go',label='g', markersize=5)
        ax1.axhline(y=1.,color='k',ls='dashed')
        ax1.axhline(y=1.1,color='k',ls='dashed')
        ax1.legend(loc='upper left')
        ax1.yaxis.set_major_locator(my_locator)
        #if myband == 'r':
        ax2.set_ylim([0.9,1.3])
        ax2.plot(normexps_r,normcounts_r,'ro',label='r', markersize=5)
        ax2.axhline(y=1.,color='k',ls='dashed')
        ax2.axhline(y=1.1,color='k',ls='dashed')
        ax2.legend(loc='upper left')
        ax2.yaxis.set_major_locator(my_locator)
        #if myband == 'i':
        ax3.set_ylim([0.9,1.3])
        ax3.plot(normexps_i,normcounts_i,'yo',label='i', markersize=5)
        ax3.axhline(y=1.,color='k',ls='dashed')
        ax3.axhline(y=1.1,color='k',ls='dashed')
        ax3.set_ylabel('Normalized counts to Exp %s 20140817' % firstExp)
        ax3.legend(loc='upper left')
        ax3.yaxis.set_major_locator(my_locator)
        #if myband == 'z':
        ax4.set_ylim([0.9,1.3])
        ax4.plot(normexps_z,normcounts_z,'mo',label='z', markersize=5)
        ax4.axhline(y=1.,color='k',ls='dashed')
        ax4.axhline(y=1.1,color='k',ls='dashed')
        ax4.legend(loc='upper left')
        ax4.yaxis.set_major_locator(my_locator)
        #if myband == 'Y':
        ax5.set_ylim([0.9,1.3])
        ax5.plot(normexps_y,normcounts_y,'co',label='Y', markersize=5)
        ax5.axhline(y=1.,color='k',ls='dashed')
        ax5.axhline(y=1.1,color='k',ls='dashed')
        ax5.legend(loc='upper left')
        ax5.yaxis.set_major_locator(my_locator)
        
        f.subplots_adjust(hspace=0)
        #plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    
        #horizontal line
        #plt.axhline(y=1.,color='k',ls='dashed')
        #plt.axhline(y=1.1,color='k',ls='dashed')
        #plt.legend(loc='upper left')
        mylabel = 'Normalized counts to Exp %s from 20140817' % firstExp
        #ax.set_ylabel(mylabel)
        myxlabel = 'Exposure Number - %s ' % firstExp    
        plt.xlabel(myxlabel)    
        #plt.title('Flat Field Count Level vs Exposures')
        outfilename = 'ExposureVsNormCounts' + '_' + ccdnum +'.pdf'
        #plt.show()
        plt.savefig(outfilename)
        plt.close()
        
        

        
    #plt.show()
    
    
    """
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
    #plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=20))

    #exp times for filters: r,z,Y = 10sec
    #exp times for filters: i = 22sec
    #exp times for filters: u,g = 30sec
    meanCountsPerSec = []
    meanFloats = []
    if infile.endswith('z.txt'):
        #theTitle = 'Counts in 10s z dome flats'
        mean = [float(x) for x in mean]
        std = [float(x) for x in std]
        #meanFloats = mean.astype(np.float)
        meanCountsPerSec = [val/10. for val in mean]
        color = 'oy'
        label = 'z'
    if infile.endswith('r.txt'):
        #theTitle = 'Counts in 10s r dome flats'
        mean = [float(x) for x in mean]
        std = [float(x) for x in std]
        #meanFloats = mean.astype(np.float)
        meanCountsPerSec = [val/10. for val in mean]
        color = 'or'
        label = 'r'
    if infile.endswith('i.txt'):
        #theTitle = 'Counts in 22s i dome flats'
        mean = [float(x) for x in mean]
        std = [float(x) for x in std]
        #meanFloats = mean.astype(np.float)
        meanCountsPerSec = [val/22. for val in mean]
        color = 'og'
        label = 'i'
    if infile.endswith('Y.txt'):
        #theTitle = 'Counts in 10s Y dome flats'
        mean = [float(x) for x in mean]
        std = [float(x) for x in std]
        #meanFloats = mean.astype(np.float)
        meanCountsPerSec = [val/10. for val in mean]
        color = 'om'
        label = 'Y'
    if infile.endswith('u.txt'):
        #theTitle = 'Counts in 30s u dome flats'
        mean = [float(x) for x in mean]
        std = [float(x) for x in std]
        #meanFloats = mean.astype(np.float)
        meanCountsPerSec = [val/30. for val in mean]
        color = 'ok'
        label = 'u'
    if infile.endswith('g.txt'):
        #theTitle = 'Counts in 30s g dome flats'
        mean = [float(x) for x in mean]
        std = [float(x) for x in std]
        #meanFloats = mean.astype(np.float)
        meanCountsPerSec = [val/30. for val in mean]
        color= 'ob'
        label = 'g'

    legend = 'Filter'
    #output pdf file name
    #tmp = infile.split('.')
    #outfilename = tmp[0] + '.pdf'
    
    #print meanCountsPerSec, type(meanCountsPerSec)
    plt.plot(xdates,meanCountsPerSec,color,label=label)
    theTitle = 'Median counts/s ugrizY'
    theYlabel = 'Median counts/sec of exposure (ADU)'
    plt.legend(loc='upper right',frameon=True,scatterpoints=1,ncol=2)
    plt.ylabel(theYlabel)
    plt.title(theTitle)
    plt.xlim(xrange)
    plt.ylim([100,3000])
    plt.gcf().autofmt_xdate()
    
    outfilename = 'mean_counts_ugrizY.pdf'
    plt.grid()
    plt.savefig(outfilename)
    print "plot Saved, waiting to display.....\n"

    plt.show()

    """
    """
    ax3 = fig.add_subplot(2,2,3,sharex=ax1)
    ax3.plot(date_to_plot,(float(calculated_data[i]['GAIN A'])/float(header_data[i]['gaina'])),'ro')
    pl.ylim=([0.7,1.3])
    
    ax4 = fig.add_subplot(2,2,4, sharey=ax3,sharex=ax1)
    ax4.plot(date_to_plot,(float(calculated_data[i]['GAIN B'])/float(header_data[i]['gainb'])),'ro')
    pl.ylim=([0.7,1.3])
    
    for ax in ax1,ax2,ax3,ax4:
        ax.grid(True)

    for l in range(61):
        if l+1 == 61:
            continue
        print "l+1", l+1
        if l <= 9:
            ax = pl.Subplot(fig,grid_spec[l,0])
            date_to_plot = datetime.datetime.strptime(str(int(calculated_data[l+1]['NIGHT'])),'%Y%m%d')
            ax.plot_date(date_to_plot,calculated_data[l+1]['GAIN A'],'ro')
            ax.plot_date(date_to_plot,header_data[l+1]['gaina'],'bo')
            if l == 0:
                pl.setp(ax.get_xticklabels(), visible=False)
            ax.grid(True)
            fig.add_subplot(ax)
        elif (l > 9) and (l <= 19):
            print "second",l-10
            second = l-10
            ax = pl.Subplot(fig,grid_spec[second,1])
            date_to_plot = datetime.datetime.strptime(str(int(calculated_data[l+1]['NIGHT'])),'%Y%m%d')
            print date_to_plot
            ax.plot_date(date_to_plot,calculated_data[l+1]['GAIN A'],'ro')
            ax.plot_date(date_to_plot,header_data[l+1]['gaina'],'bo')
            if l ==10:
                pl.setp(ax.get_xticklabels(), visible=False)
            ax.grid(True)
            fig.add_subplot(ax)
        elif (l > 19) and (l <= 29):
            print "third",l-20
            third = l-20
            ax = pl.Subplot(fig,grid_spec[third,2])
            date_to_plot = datetime.datetime.strptime(str(int(calculated_data[l+1]['NIGHT'])),'%Y%m%d')
            print date_to_plot
            ax.plot_date(date_to_plot,calculated_data[l+1]['GAIN A'],'ro',)
            ax.plot_date(date_to_plot,header_data[l+1]['gaina'],'bo')
            if l == 19:
                pl.setp(ax.get_xticklabels(), visible=False)
            ax.grid(True)
            fig.add_subplot(ax)
            
            
        #print calculated_data
        #print header_data
        #pl.plot(xdates[i],gain_A)
        #pl.plot(header_ccd,header_gain_A)
        #pl.ylim(ymax=6.0,ymin=0)
        #pl.xlim(xmin=1,xmax=62)
    pl.show()
    """
if  __name__ == '__main__':
    main()
    