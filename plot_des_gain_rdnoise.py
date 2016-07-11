#!/usr/bin/env python
"""Plot the gain and rdnoise for each ccd and compares with
the values in the image headers.

.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>

  Needs to setup the following before running:
  matplotlib
  
"""
# imports
import coreutils.desdbi
import argparse
import os
import csv
from matplotlib import pyplot as pl
from matplotlib.gridspec import GridSpec
import datetime
import matplotlib.dates as mdates


parser = argparse.ArgumentParser(description="Plot gain and rdnoise for each ccd for every night there is a measurement")
parser.add_argument("-l", "--list", required=True, dest="list_of_files", help="List of files for each night with the gain and rdnoise for each ccd, comma-separated")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")

args = parser.parse_args()
all_files = args.list_of_files.split(',')
verbose = args.verbose

def main():

    
    dates = []
    #find all the dates available from the file names:
    for files in all_files:
        tmp = files.split('.')
        dates.append(tmp[0].split('_')[2])

    xdates = [datetime.datetime.strptime(str(int(date)),'%Y%m%d') for date in dates]
    
    #fig_gains = pl.figure(figsize=(8, 8))
    fig_gains, (ax1, ax2) = pl.subplots(1,2, sharex=True, sharey=True)

    #fig_gains_ratios = pl.figure(figsize=(8, 8))
    #grid_spec = GridSpec(10, 3, wspace=0.0, hspace=0.0)

       
    all_dates = []
    all_calculated_data = []
    all_head_data = []
    for i,name in enumerate(all_files):
        tmp = name.split('.')
        date = tmp[0].split('_')[2]

        
        #gain and rdnoise from image header
        file_data_header = 'gain_rdnoise_db_image_'+ date + '.csv'
        
        with open(name,'rb') as data_file:
            for i in range(4):
                next(data_file)
            #reader = csv.reader(data_file)
            reader = csv.DictReader(data_file, delimiter=',')
            calculated_data = [r for r in reader]

        #read the data from file with header values
        with open(file_data_header,'rb') as values_header:
            reader = csv.DictReader(values_header, delimiter=',')
            header_data = [r for r in reader]
 
        #this can be used to determine the range of dates
        all_dates.append(date)
        all_calculated_data.append(calculated_data)
        all_head_data.append(header_data)
        
    print "all_dates", all_dates
    #all_dates.append('20121101')
    #all_dates.append('20130205')
    cal_gain_A_1 = []
    cal_gain_B_1 = []
    cal_rdnoise_A_1 = []
    cal_rdnoise_B_1 = []

    cal_gain_A_31 = []
    cal_gain_B_31 = []
    cal_rdnoise_A_31 = []
    cal_rdnoise_B_31 = []

    cal_gain_A_45 = []
    cal_gain_B_45 = []
    cal_rdnoise_A_45 = []
    cal_rdnoise_B_45 = []

    head_gain_A_1 = []
    head_gain_B_1 = []
    head_rdnoise_A_1 = []
    head_rdnoise_B_1 = []

    head_gain_A_31 = []
    head_gain_B_31 = []
    head_rdnoise_A_31 = []
    head_rdnoise_B_31 = []

    head_gain_A_45 = []
    head_gain_B_45 = []
    head_rdnoise_A_45 = []
    head_rdnoise_B_45 = []

    date_to_plot = []
    #for date in all_dates:
        #print date
    #date_to_plot = [datetime.datetime.strptime(d,'%Y%m%d').date() for d in all_dates]
        #date_to_plot.append(datetime.datetime.strptime(str(int(date)),'%Y%m%d'))
        
    #print date_to_plot
    #print all_calculated_data[0]
    #print len(all_calculated_data[0])
    #print all_calculated_data[1]
    #print len(all_calculated_data[1])
    
    for item in all_calculated_data:
        #print all_calculated_data[i]
        cal_gain_A_1.append(item[0]['GAIN A'])
        cal_gain_A_31.append(item[30]['GAIN A'])
        cal_gain_A_45.append(item[44]['GAIN A'])
        cal_gain_B_1.append(item[0]['GAIN B'])
        cal_gain_B_31.append(item[30]['GAIN B'])
        cal_gain_B_45.append(item[44]['GAIN B'])

    for item in all_head_data:
        head_gain_A_1.append(item[0]['gaina'])
        head_gain_A_31.append(item[30]['gaina'])
        head_gain_A_45.append(item[44]['gaina'])
        head_gain_B_1.append(item[0]['gainb'])
        head_gain_B_31.append(item[30]['gainb'])
        head_gain_B_45.append(item[44]['gainb'])
    
    print all_dates
    print head_gain_A_31
    print head_gain_B_31
    print cal_gain_A_31
    print cal_gain_B_31
    
    print "GAIN A, HEAD A, GAIN B, HEAD B,ratio A (cal/head), ratio B (cal/head)"
    for i in range(len(all_dates)):
        
        print "CCD 1 date: %s %4.2f  %4.2f  %4.2f  %4.2f %4.2f %4.2f" % (all_dates[i],float(cal_gain_A_1[i]),float(head_gain_A_1[i]),float(cal_gain_B_1[i]),float(head_gain_B_1[i]), float(cal_gain_A_1[i])/float(head_gain_A_1[i]),float(cal_gain_B_1[i])/float(head_gain_B_1[i]) )

    for i in range(len(all_dates)):

        print "CCD 31 date: %s %4.2f  %4.2f  %4.2f  %4.2f %4.2f %4.2f" % (all_dates[i],float(cal_gain_A_31[i]),float(head_gain_A_31[i]),float(cal_gain_B_31[i]),float(head_gain_B_31[i]),float(cal_gain_A_31[i])/float(head_gain_A_31[i]),float(cal_gain_B_31[i])/float(head_gain_B_31[i]) )

    
        #print "date: %s gain A   ccd31  %4.2f  gain B   %4.2f" % (all_dates[i],float(cal_gain_A_31[i]),float(cal_gain_B_31[i]))
    
    for i in range(2):
        
        
        # use a more precise date string for the x axis locations in the
        # toolbar
    
        #print calculated_data[i]['NIGHT']
        date_min= datetime.datetime.strptime('20121101','%Y%m%d').date()
        date_max= datetime.datetime.strptime('20130205','%Y%m%d').date()


        print len(cal_gain_A_1)

        #print "xdates",date_to_plot
        print all_dates
        
        date_to_plot = [mdates.date2num(date) for date in xdates]
        if i == 0:
            #ccd 1
            #for j in range(all_dates):
                #date_to_plot = datetime.datetime.strptime(thedates,'%Y%m%d').date()
            print date_to_plot
            
            #ax = fig_gains.add_subplot(1,1,1)
            ax1 = fig_gains.add_subplot(2,1,1)
            ax1.plot(date_to_plot,cal_gain_A_1,'ro')
            ax1.plot(date_to_plot,head_gain_A_1,'bo')
            ax1.set_xlim([date_min,date_max])
            ax1.set_ylim([4.2,4.6])
            ax1.set_ylabel('Gain e/ADU')
            fig_gains.suptitle('Gain Amp A and Amp B, CCD 1')
            
            ax2 = fig_gains.add_subplot(2,2,1,sharey=ax1)
            ax2.plot(date_to_plot,cal_gain_B_1,'ro')
            ax2.plot(date_to_plot,head_gain_B_1,'bo')
                
            ax1.grid(True)
            ax2.grid(True)
            fig_gains.fmt_xdata = mdates.DateFormatter('%Y%m%d')
            fig_gains.autofmt_xdate(bottom=0.2,rotation=45,ha='right')
        #pl.setp(ax1.get_xticklabels(), visible=True)
            #ax.set_ylabel('Gain e/ADU')
            #ax.set_title('Gain CCD1')
        
        print i
    #fig_gains.fmt_xdata = mdates.DateFormatter('%Y%m%d')
    #fig_gains.autofmt_xdate()
    #pl.savefig('name.jpg')
    pl.show()

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
