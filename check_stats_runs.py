#!/usr/bin/env python
"""
Check the status of runs being processed suing desstat

.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>

This code check how many runs are currently running


"""

# imports
import sys
import os
import subprocess
import datetime,time
    
#    Turn on svn keyword substitution for Revision (above) with:
#    bash$ svn propset svn:keywords "Revision" foo_quality.py
__version__ = "$Revision: 187 $"

               


def main():


    try:
        des_home = os.environ["des_home"]
    except KeyError:
        des_home  = None

    if des_home is None:
        print "setting up desdmsoft current"
        #set_deshome = subprocess.Popen('$EUPS_DIR/bin/eups_setup desdmsoft 7.2.9+1', shell=True)
        os.system('./setup_desdmsoft.csh')

    p = subprocess.Popen('desstat', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    for line in p.stdout.readlines():
        print line

    start = str(datetime.datetime.now())
    now = str(datetime.datetime.now())
    end = '2013-09-11 07:00:00.000000'
    
    current_time = now.split(' ')
    current_hrs = current_time[1].split(':')[0]
    print current_hrs
    
    
    while int(current_hrs) <= 23:
        time.sleep(20)
        p = subprocess.Popen('desstat', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        with open('destats_check.tmp','a') as outfile:    
            unix_time = str(time.time())
            print unix_time
            mytime = str(datetime.datetime.now())
            outfile.write(mytime+'\n')
            outfile.write(unix_time)
            
            for line in p.stdout.readlines():
                #if line[-1] == '\n':
                #    line = line[:-1]
                if line[0] == '\n':
                    line = line[0:]
                outfile.write(line)
    #outfile.close()
    
    """
    #check the tmp file to determine if we are ready to submit
    with open('destats_check.tmp','r') as input:
        useful_lines = input.readlines()[3:]
        for line in useful_lines:
            #parse all the lines:
            cols = line.split()
            print 'cols',cols
        
    """ 
    """
    What I need to know is how many machines firstcut takes on average at each block
    """
        

if  __name__ == '__main__':
    main()


