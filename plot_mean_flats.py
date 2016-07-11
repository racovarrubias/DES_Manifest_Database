#!/usr/bin/env python
"""
Plots the mean flat counts/s for exposiure per night.
It uses the output file from des_gain_rdnoise.py

.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>

  Needs to setup the following before running:
  matplotlib
  
"""
# imports
#import coreutils.desdbi
import argparse
import os
import csv
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import datetime as dt
import matplotlib.dates as mdates
import numpy as np


parser = argparse.ArgumentParser(description="Plot gain and rdnoise for each ccd for every night there is a measurement")
parser.add_argument("-f", "--file", required=True, dest="infile", help="Input file with mean counts and std per night (output from des_gain_rdnoise.py)")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")

args = parser.parse_args()
files = args.infile.split(',')
verbose = args.verbose

def main():


    # Generate a series of dates (these are in matplotlib's internal date format)
    dates = mdates.drange(dt.datetime(2013, 01, 01), dt.datetime(2013,12,10),dt.timedelta(weeks=6))
    
    #find all the dates available from the file names:
    for infile in files:
        dates = []
        mean = []
        std = []
        with open(infile,'rb') as inputfile:
        
            for line in inputfile:
                #line = inputfile.readline()
                line = line.rstrip()
                values = line.split(' ')
                mean.append(values[0])
                std.append(values[1])
                dates.append(values[3])
        

        xdates = [dt.datetime.strptime(str(int(date)),'%Y%m%d') for date in dates]
        xrange = [dt.datetime.strptime(str(int(20130901)),'%Y%m%d'),dt.datetime.strptime(str(int(20131210)),'%Y%m%d')]
        #print xdates
        #print dates

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
