#!/usr/bin/env python
"""
Looks for json manifest file sent via DTS, reads it and then 
looks in the exposure table to check if all the fileds for the field
observed arrived.
Submits firstcut_snse


.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>


This code submit firstcut_snse.des for every field. It needs to find the 
exposure.id for each of the exposures int the OBS_SET_EXPOSURE that 
it's going to process.
The tables that contain the information about the exposures taken during
a night are in the following tables.

   Todd created these tables in destest DB:
    TABLE_NAME
    ------------------------------
    OBS_SET
    =======
    ID
    TYPE
    NITE
    CREATED_DATE
    

    OBS_SET_TYPE
    ============
    TYPE
    
    This table contains all the exposures observed for a given OBS_SET_ID 
    OBS_SET_EXPOSURE
    ================
    OBS_SET_ID
    EXPOSURENAME
    

OBS_SET   /* an unordered group of exposures taken during an observing night */
-----------------------------
ID  INTEGER  NOT NULL    /* auto increment ID from sequence */
TYPE  VARCHAR2(16)   /* from a set of valid types defined in another table TBD */
NITE   VARCHAR2(8)   /* the NITE of observing, in YYYYMMDD form  */
CREATED_DATE  DATE  /* auto-filled date and time the row was created  */


OBS_SET_EXPOSURE  /* identifies the OBS_SET to which an EXPOSURE belongs */
--------------------------
OBS_SET_ID  INTEGER  NOT NULL  /* the OBS_SET  */
EXPOSURENAME  VARCHAR2(100)  NOT NULL  /* the name of the exposure file */

OBS_TYPE /* SN-C1, SN-C2, SN-C3-g, SN-C3-i, SN-C3-r, SN-C3-z, SN-E1, SN-E2, SN-S1 ...
 

MANIFEST_SUMMARY  /* a view that aggregates sets and exposures from a NITE  */
----------------------------
NITE
NUM_OBS_SETS
NUM_EXPOSURES
NUM_INGESTED


The code reads a template file stored in submitfiles_templates directory
and it modifies the critical variables for each run. 
the variables to modify in the SNeSE pipeline are:

snese.des

PIPELINE = snese
JIRA_ID = DESOPS-400
EVENT_TAG = snese_20121216_attemp1
RUN_COMMENT = snese_20121216_attemp1
NITE = 20121216
PROJECT = OPS

After SNeSE is finished, it needs to start the diffim software.

diffim submit file:

QUERY_RUN = 20130115124130_20130114   # This run is the run from firstcut_snse

NITE = 20130114
RUN_COMMENT = Diffimg_reprocessing_2013_20130114_E1_cmp2_attempt1
EVENT_TAG = diffimg_20130114_E1_cmp2_fermi_attemp1

SNPOINTQ = "SN-E1"
JIRA_ID = DESOPS-403

"""

# imports
import coreutils.desdbi
import argparse
import os
import sys
import time
from datetime import datetime, date, time, timedelta
import string
import csv

parser = argparse.ArgumentParser(description="Check the SN Manifest json file sent for every SNe block and when all files arrive submit the firstcut_snse job")
parser.add_argument("-n", "--nite", required=True, dest="night", help="Night to process")
parser.add_argument("-p", "--project", required=False, dest="project", default='OPS', help="Project to process data. Default is OPS")
parser.add_argument("-l", "--list_type", required=False,action="store_true", help="Observation type. Run submit_sne_run.py -list_types to list all valid types")
parser.add_argument("-j", "--jira_id", required=True, dest="jira_id", help="Optional Jira id. If not given, it's transfered from template submit file")
parser.add_argument("-e", "--event_tag", required=True, dest="event_tag", help="Optional Event. If not given, it's transfered from template submit file")
parser.add_argument("-c", "--run_comment", required=True, dest="run_comment", help="Optional run_comment. If not given, it's transfered from template submit file")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")
#parser.add_argument("-i", "--idlist", required=True, dest="idlist", help="Path to ID list of images to process")
#parser.add_argument("-t", "--obstype", required=True, dest="obstype", help="Observation type. Run submit_sne_run.py -l to list all valid types")

args = parser.parse_args()
night = args.night
#idlist = args.idlist
#obstype = args.obstype
list_type = args.list_type
project = args.project
jira_id = args.jira_id
event = args.event_tag
run_comment = args.run_comment
verbose = args.verbose

    
#    Turn on svn keyword substitution for Revision (above) with:
#    bash$ svn propset svn:keywords "Revision" foo_quality.py
__version__ = "$Revision: 187 $"




def check_if_exp_ingested(exposure, night, cur):
    
    #write a query to see if the file has already been ingested in the database.
    queryitems = ["ose.exposurename", "os.type", "e.exposurename"]
    #queryitems = ["exposurename, nite, band, field, sn_sequence_status"]
    querylist = ",".join(queryitems)

    
    queryitems = ["id,exposurename, nite"]
    querylist = ",".join(queryitems)

    
    query = """select %s from exposure
            where exposurename like '%s%%'
            and exposuretype='src'
            and project = 'DTS'
            and nite = '%s' """ % (querylist, exposure, night)
    
    print query    
    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)

    return cur

def connectDB(dbase):
    """
    Connect to database, query and extract data from table, and close connection to database
    dbase can be: 
    db-destest
    db-desoper
    """
    try:
        desdmfile = os.environ["des_service"]
    except KeyError:
        desdmfile  = None


    dbh = coreutils.desdbi.DesDbi(desdmfile,dbase)
    if dbh.is_postgres():
        print 'Connected to postgres DB'
    elif dbh.is_oracle():
        print 'Connected to oracle DB'
        print 'which_services_file = ', dbh.which_services_file()
        print 'which_services_section = ' , dbh.which_services_section()
        
    cur = dbh.cursor()
        
    return cur
               

def create_id_list(night=night):
    """
    Creates a file with all the exposure.id for all the SNe exposures to process for given night
    :param: night: Night to process
    """
    mycur = connectDB("db-desoper")

    id_query = """select id from exposure 
                    where nite = '%s'
                    and object like '%%SN-%%'
                    and exptime > '15'
                    order by exposurename""" % (night)

 
    mycur.arraysize = 1000 # get 1000 at a time when fetching
    mycur.execute(id_query)
    outfilename = 'SN_' + night + '.ids'
    rowcount = 0
    with open(outfilename,'w') as outfile:       
        for item in mycur:
            ids = item[0]
            print ids
    sys.exit(0)

    
    

def submit_firstcut_snese(night=night, jira=None, event=None, project=None, run_comment = None, idfile = None):
    """
    Reads submit for single epoch SNe processing file from submit directory and updates 
    parameters
    :param night: Night to process
    :param jira_id: Optional Jira_id where notes are taken
    :param project: Optional project for processing. Default is OPS, but can be ACT for testing
    
    """
    
    submit_file = 'firstcut_snese_'+night+'.des'
    #open the template file for firstcut_snese
    with open(submit_file,'w') as outfile:
        with open('firstcut_snese.des','r') as temp_file:
            for line in temp_file:
                if line.endswith('.ids>> '):
                    if idfile == None:
                        print "Must provide a file with exp.id values to process"
                        sys.exit(0)
                    else:
                        line.replace('ID_LIST_FILE.ids',idlist)
                        outfile.write(line)
                if line.startswith('JIRA_ID = '):
                    if jira == None:
                        print "Using default value from template file", line
                        continue
                    else:
                        line = line[0:len('JIRA_ID = ')] + jira
                        outfile.write(line)
                if line.startswith('EVENT_TAG = '):
                    if event == None:
                        print "Using default value from template file", line
                        continue
                    else:
                        line = line[0:len('EVENT_TAG = ')] + event
                        outfile.write(line)
                if line.startswith('RUN_COMMENT = '):
                    if run_comment == None:
                        print "Using default value from template file", line
                        continue
                    else:
                        line = line[0:len('RUN_COMMENT = ')] + run_comment
                        outfile.write(line)
                if line.startswith('NITE = '):
                    if night == None:
                       print "Must provide a night to process"
                       sys.exit(0)
                    else: 
                        line = line[0:len('NITE = ')] + night
                        outfile.write(line)
                if line.startswith('PROJECT = '):
                    if project == None:
                        print "Using ACT project..."
                        line = line[0:len('PROJECT = ')] + 'ACT'
                    else:
                        line = line[0:len('project = ')] + night
                        outfile.write(line)
                
                    
                
    
    
def submit_firstcut(night, jira_id):
    """
    Reads submit file for firstcut from submit directory and updates 
    parameters
    :param night: Night to process
    :param jira_id: Optional Jira_id where notes are taken
    :param project: Optional project for processing. Default is OPS, but can be ACT for testing
    
    """
    


def submit_diffim(red_run, field, jira_id):
    """
    Reads submit file for difference image from submit directory and updates 
    parameters
    :red_run : red run from firstcut_snese
    :field : field to do diffim that is present in red_run 
    :param jira_id: Optional Jira_id where notes are taken
        
    """
    
    
def check_running_runs():
    """
    Check for status of run currently running on iforge and desjobc
    """
    


def main():


    create_id_list(night)
    
    #list all observation type present in OBS_SET_TYPE table.
    if list_type:

        mycur = connectDB("db-destest")
        
        
        print "Current type available:"
        type_query = "select distinct TYPE from OBS_SET_TYPE"
        mycur.arraysize = 1000 # get 1000 at a time when fetching
        mycur.execute(type_query)
        for item in mycur:
            print item[0]
        sys.exit(0)


    
    """
    Todd created these tables in destest DB:
    TABLE_NAME
    ------------------------------
    OBS_SET:
    ID
    TYPE
    NITE
    CREATED_DATE
    

    OBS_SET_TYPE:
    TYPE
    
    This table contains all the exposures observed for a given OBS_SET_ID 
    OBS_SET_ID
    EXPOSURENAME
    
    OBS_SET_TYPE: This table contains the type of observations done at each night. The types
    can be: survey, dflat-g, dflat-r, ..., SN-X3-g,SN-X3-r, .., SN-C1,SN-C2...

    To select the exposure names for a given type you can do the following query:
    select ose.exposurename,os.type from obs_set_exposure ose, obs_set os where os.type='survey' and ose.obs_set_id = os.id and os.nite='20121017';
    
    This query will tell you if the exposures have arrived and ingested into the DB at NCSA:
    select ose.exposurename,os.type,e.exposurename from exposure e,obs_set_exposure ose, obs_set os where os.type='survey' and ose.obs_set_id = os.id and os.nite='20121017' and e.exposurename = ose.exposurename;
    
    A summary of all observations can be seen in the view table:
    MANIFEST_SUMMARY
    NITE  
    NUM_OBS_SETS 
    NUM_EXPOSURES 
    NUM_INGESTED  

    """

    #try:
    #    desdmfile = os.environ["des_services"]
    #except KeyError:
    #    desdmfile = None
    #dbh = coreutils.desdbi.DesDbi(desdmfile,"db-destest")
    #cur = dbh.cursor()
    cur = connectDB("db-destest")
    
    queryitems = ["ose.exposurename", "os.type", "e.exposurename"]
    #queryitems = ["exposurename, nite, band, field, sn_sequence_status"]
    querylist = ",".join(queryitems)

    #check current date to determine what date we want to process
    #The data to process has to be from dates within one day of current date.
    #We don't want to process nights from more than 1 day, assuming we have 
    #been process nights continuosly.
    if night is None:
        #This will execute a query on the SN_MANIFEST table for data within the day the script is execute and 
        #and the previus day.
        today = datetime.now()
        yesterday = today - timedelta(days=1)
        
        today_year = today.year
        today_month = today.month
        today_day = today.day
        
        yesterday_year = yesterday.year
        yesterday_month = yesterday.month
        yesterday_day = yesterday.day
        
        
        if today_month < 10:
            str_today_month = '0'+str(today_month)
        else:
            str_today_month = str(today_month)
        if today_day < 10:
            str_today_day = '0'+str(today_day)
        else:
            str_today_day = str(today_day)
    
        if yesterday_month < 10:
            str_yesterday_month = '0'+str(yesterday_month)
        else:
            str_yesterday_month = str(yesterday_month)
        if yesterday_day < 10:
            str_yesterday_day = '0'+str(yesterday_day)
        else:
            str_yesterday_day = str(yesterday_day)
            
        #write yesterday date as string so I can compare with the database.
        str_yesterday = str(yesterday_year) + str_yesterday_month + str_yesterday_day
        str_today = str(today_year) + str_today_month + str_today_day
        
        print str_today, str_yesterday

        query = """select %s 
        from exposure e,obs_set_exposure ose, obs_set os 
        where os.type = '%s'
        and ose.obs_set_id = os.id 
        and os.nite > '%s'
        and os.nite < '%s'
        and e.exposurename = ose.exposurename
        order by e.exposurename""" % (querylist,  str_yesterday, str_today)
        
        print query
    
    else:
        query = """select %s 
        from exposure e,obs_set_exposure ose, obs_set os 
        where os.type='survey' 
        and ose.obs_set_id = os.id 
        and os.nite = '%s'
        and e.exposurename = ose.exposurename
        order by e.exposurename""" % (querylist, night)
    
        print query

    #execute the query
    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)
    
    
    #Preparing the files to submit
    obs_exposures = [] #exposure taken at the mountain
    ingest_exposures = [] #exposures that have arrived and ingested at NCSA
    obs_type = [] #exposure type
    for item in cur:
        #fileds in the query: ose.exposurename", "os.type", "e.exposurename"
        obs_exposures.append(item[0])
        ingest_exposures.append(item[2])
        obs_type.append(item[1])
        
    #print obs_exposures
    #print ingest_exposures
    #print obs_type
    
                
    #let's check if the exposure has arrived and ingested
    #result = check_if_exp_ingested(exposure, obs_night, mycur)
    #print result
    #submit_firstcut_snese(night, jira=jira_id, event=event_tag, project=project, run_comment = comment)

if  __name__ == '__main__':
    main()


