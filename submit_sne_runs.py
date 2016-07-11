#!/usr/bin/env python

"""
Submits snese and Diffim runs


.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>


This code submit snese.des for every field. It needs to find the 
exposure.id for each of the exposures ingested in the OBS_SET_EXPOSURE, table that
is filled with the ingestManifest.py file.
The tables that contain the information about the exposures taken during
a night are:


OBS_SET   /* an unordered group of exposures taken during an observing night */
-----------------------------
ID  INTEGER  NOT NULL    /* auto increment ID from sequence */
TYPE  VARCHAR2(16)   /* from a set of valid types defined in another table TBD */
NITE   VARCHAR2(8)   /* the NITE of observing, in YYYYMMDD form  */
CREATED_DATE  DATE  /* auto-filled date and time the row was created  */


OBS_SET_EXPOSURE  /* identifies the OBS_SET to which an EXPOSURE belongs */
--------------------------
OBS_SET_ID  INTEGER  NOT NULL  /* the OBS_SET  */
EXPOSURENAME  VARCHAR2(100)  NOT NULL  /* the name of the exposure file */   This relates the exposurename (DECam_00123456.fits) to the obs_type


OBS_SET_TYPE /* Table with a defined set of exposure types
---------------------------
TYPE    VARCHAR2(100)        Observation type (dflat-r, SN-C1, Zero, ....)


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


The code reads a template file stored in submitfiles_templates directory
and it modifies the critical variables for each run. 
the variables to modify in the snese, and diffim
pipelines are:

snese.des

PIPELINE = snese
JIRA_ID = DESOPS-400
EVENT_TAG = snese_20121216_attemp1
RUN_COMMENT = snese_20121216_attemp1
NITE = 20121216
PROJECT = OPS

After snese is finished, it needs to start the diffim software.

diffim submit file key values to change:

QUERY_RUN = 20130115124130_20130114   # This run is the run from snese

NITE = 20130114
RUN_COMMENT = Diffimg_20131109_S1_attemp1
EVENT_TAG = Diffimg_20131109_S1_attemp1

SNPOINTQ = "SN-E1"
JIRA_ID = DESOPS-403

Examples:

submit_sne_runs.py -n 20131226 -p OPS -t SN-S1 -v -et None

"""

# imports
import coreutils.desdbi
import argparse
import os
import sys
import time
from datetime import datetime, date, time, timedelta
import string
import json
import subprocess

parser = argparse.ArgumentParser(description="Check the SN Manifest json file sent for every SNe block and when all files arrive submit the firstcut_snse job")
parser.add_argument("-n", "--nite", required=True, dest="night", help="Night to process")
parser.add_argument("-p", "--project", required=True, dest="project", default='OPS', help="Project to process data. Default is OPS")
parser.add_argument("-i", "--idlist", required=False, dest="idlist", help="Path to ID list of images to process")
parser.add_argument("-t", "--obstype", required=True, dest="obstype", help="Observation type. Run submit_sne_run.py -l to list all valid types")
parser.add_argument("-l", "--list_type", required=False,action="store_true", help="List all observation type. Run submit_sne_run.py -list_types to list all valid types")
parser.add_argument("-bc", "--biascor", required=False, dest="biascor", help="Optional biascor run. If not given, it's transfered from template submit file")
parser.add_argument("-fc", "--flatcor", required=False, dest="flatcor", help="Optional flatcor run. If not given, it's transfered from template submit file")
parser.add_argument("-jid", "--jira_id", required=False, dest="jira_id", help="Optional Jira id. If not given, it's transfered from template submit file")
parser.add_argument("-et", "--eventtag", required=False, dest="event", help="Optional Event. If not given, it's transfered from template submit file")
parser.add_argument("-rc", "--run_comment", required=False, dest="run_comment", help="Optional run_comment. If not given, it's transfered from template submit file")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")

args = parser.parse_args()
night = args.night
idlist = args.idlist
obstype = args.obstype
list_type = args.list_type
biascor = args.biascor
flatcor = args.flatcor
project = args.project
jira_id = args.jira_id
event = args.event
run_comment = args.run_comment
verbose = args.verbose

    
#    Turn on svn keyword substitution for Revision (above) with:
#    bash$ svn propset svn:keywords "Revision" foo_quality.py
__version__ = "$Revision: 187 $"


def list_observations_types(mycur):
    """ 
    List all the observation types from DB table OBS_SET_TYPE
    :param mycur: Database connection
    
    """
        
    print "Current type available:"
    type_query = "select distinct TYPE from OBS_SET_TYPE"
    mycur.arraysize = 1000 # get 1000 at a time when fetching
    mycur.execute(type_query)
    for item in mycur:
        print item[0]

    sys.exit(0)
    

def check_observed_field_json(json_file):
    """
    Looks for the manifest file sent via DTS. this file contains information
    about the exposures taken for a given SN field
    """
    
    #open json file with exposures information
    #I am interested in just the expid. I need to check if the expid's were ingested.
    #I also need to select the id from the exposure table corresponding to this expid's
    with open(json_file, 'rb') as jf:
        data = json.load(jf)  #data is a dictionary
        
 
 
   

    
def check_exposures_in_manifest(night,querylist,cur):
    
    
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
        order by e.exposurename""" % (querylist,  obstype, str_yesterday, str_today)
        
        print query
    
    else:
        query = """select %s 
        from exposure e,obs_set_exposure ose, obs_set os 
        where os.type='%s' 
        and ose.obs_set_id = os.id 
        and os.nite = '%s'
        and e.exposurename = ose.exposurename
        order by e.exposurename""" % (querylist, obstype, night)
    
        print query

        #execute the query
    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)

        #Preparing the files to submit
    exposures_in_manifest = [] #exposure in manifest file (exposure that was sent from CTIO)
    ingested_exposures = [] #exposures that have arrived and ingested at NCSA
    obs_type = [] #exposure type
    for item in cur:
        exposures_in_manifest.append(item[0])
        ingested_exposures.append(item[2])
        obs_type.append(item[1])
    
    if len(exposures_in_manifest) == 0:
        print "No exposures found for %s field for night %s" % (obstype, night)
        sys.exit(0)
        
    return exposures_in_manifest, ingested_exposures, obs_type


def manifest_exposures_equal_to_ingested_exposures(exposure_manifest,ingested_exposures,verbose=None):
    """
    Check if all exposures in manifest files have been ingested in the DB

    :param exposure_manifest:  List with exposurenames in manifest
    :param ingested_exposures: List with all exposurename found in DB for current type in manifest file
    """

    if set(exposure_manifest) == set(ingested_exposures) and len(exposure_manifest) == len(ingested_exposures):
        all_exposures_are_ingest = True
    else:
        all_exposures_are_ingest = False
        
    if len(exposure_manifest) < len(ingested_exposures):
        print "The number of ingested exposures if larger than the number of exposures in Manifest file"
        print "Exiting code...."
        sys.exit(0)
    
    if verbose:
        print "Exposures in Manifest File: \n", exposure_manifest
        print "Exposures Ingested: \n", ingested_exposures
        
    return all_exposures_are_ingest
    
    
    
def create_idlist_for_submit_file(exposures, night,cur,exposure_type,verbose=None):
    """
    Check if the exposure has been ingested and creates the idlist with
    all exposures id, expect the pointing exposure (usually the first 10 sec
    exposure in the manifest file.
    
    :param exposure: exposurename from exposure table
    :param night: night from where to get the exposurename. This is defined as the night when observations began (local night time)
    :param exposure_type: exposure type (SN-C1, SN-X3-g, ...)
    :param cur: database cursor
    :param verbose: optional verbose.
    
    """
    root_idlist_path = '/home/ricardoc/DESDM/Y1-SNe/snese/ids/' 
    
    #write a query to see if the file has already been ingested in the database.
    #queryitems = ["ose.exposurename", "os.type", "e.exposurename"]
    #queryitems = ["exposurename, nite, band, field, sn_sequence_status"]
    #querylist = ",".join(queryitems)

    
    queryitems = ["id, exptime, exposurename, nite"]
    querylist = ",".join(queryitems)


    #create idlist
    idlist_filename = exposure_type + '_' + night + '.ids'
    full_path_filename = os.path.join(root_idlist_path,idlist_filename)
    if verbose:
        print "id filename is: %s \n" % full_path_filename
        
    if not os.path.exists(full_path_filename):
        with open(full_path_filename,'w') as outfile:
            outfile.write("id=")
    
            for i,each_exposure in enumerate(exposures):
                query = """select %s from exposure
                        where exposurename like '%s%%'
                        and exposuretype='src'
                        and project = 'DTS'
                        and nite = '%s' """ % (querylist, each_exposure, night)
    
                if verbose:
                    print query    
        
                cur.arraysize = 1000 # get 1000 at a time when fetching
                cur.execute(query)
    
                for item in cur:
                    #pointing exposure doesn't have to be processed, so it is excluded from idlist.
                    if item[1] == 10.0:
                        continue
                    print item
                    if i == len(exposures) - 1:
                        outfile.write("%s" % item[0])
                    else:
                        outfile.write("%s," % item[0])

    else:
        print "File already exists. Continue..."           
        
         
    success_creating_list = True
    
    #check if the exposure is in DTS directory and ingested in DB
    #for item in cur:
    #    if verbose:
    #        print "exposure is in DTS directory:  ", item

    return success_creating_list, full_path_filename

def connectDB(dbase):
    """
    Connect to database, query and extract data from table, and close connection to database
    dbase can be: 
    db-destest
    db-desoper
    
    :param dbase: Database name to query. db-destest or db-desoper
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
        print 'which_services_section = ', dbh.which_services_section()
        
    cur = dbh.cursor()
        
    return cur
               

def create_submit_snese(night=night, jira=None, event=None, project=None, run_comment = None, idfile = idlist, type = None, biascor = None, flatcor = None):
    """
    Using a template submit file and creates a snese submit file for  given night.
    Reads submit for single epoch SNe processing file from submit directory and updates 
    parameters
    :param night: Night to process
    :param jira_id: Optional Jira_id where notes are taken
    :param event: EVENT_TAG keyword with comment to identify run processing
    :param run_comment: RUN_COMMENT keyword with comment to identify run processing
    :param project: Optional project for processing. Default is OPS, but can be ACT for testing
    :param idfile: File with exposureid from ingested files into the DB. These files will be processed
    :param type: Observation type of exposures. 
    :param biascor: biascor run to use for processing exposures
    :param flatcor: flatcor run to use for processing expoures
    """
    
    
    sneseTemplateFilePath = '/home/ricardoc/DESDM/Y1-SNe/snese/'

    #Don't apply fringe correction to bands: g,r,i
    #Apply fringe correction to bands: z
    """
    If field is C1,C2,X1,X2,S1,S2,E1,E2
    then use FullblockList
    if field is X3-g,X3-r, X3-i, C3-g, C3-r, C3-i
    then use blockListNoFringeCor
    if field is X3-z, C3-z 
    then use blockListFringeCor
    """
    
    if ('SN-C1' in idfile) or ('SN-C2' in idfile) or ('SN-S1' in idfile) or ('SN-S2' in idfile) or ('SN-E1' in idfile) or ('SN-E2' in idfile):
        blockList = 'crosstalk,imcorrect_nofringecor,imcorrect_fringecor,astrorefine,ingest_scampqa,make_bleedmask,make_craymask,compress_files'
    elif ('X3-g' in idfile) or ('X3-r' in idfile) or ('X3-i' in idfile) or ('C3-g' in idfile) or ('C3-r' in idfile) or ('C3-i' in idfile):
        blockList = 'crosstalk,imcorrect_nofringecor,astrorefine,ingest_scampqa,make_bleedmask,make_craymask,compress_files'
    elif ('X3-z' in idfile) or ('C3-z' in idfile):
        blockList = 'crosstalk,imcorrect_fringecor,astrorefine,ingest_scampqa,make_bleedmask,make_craymask,compress_files'
    
    print "block list to use is  ", blockList
    
    submit_file = os.path.join(sneseTemplateFilePath,'snese_'+ type + '_' + night + '.des')
    template_submit_file = os.path.join(sneseTemplateFilePath,'snese_Y1_template.des')
    #open the template file for firstcut_snese
    with open(submit_file,'w') as outfile:
        with open(template_submit_file,'r') as temp_file:
            for line in temp_file:
                #line = line.rstrip()
                if line.endswith('ID_LIST_FILE.ids>>\n'):
                    if idfile == None:
                        print "Must provide a file with exp.id values to process"
                        sys.exit(0)
                    else:
                        line = line[0:len('<<include ')] + idfile + '>>\n'
                        outfile.write(line)
                elif line.startswith('block_list = '):
                    line = line[0:len('block_list = ')] + blockList + '\n'
                    outfile.write(line)
                elif line.startswith('JIRA_ID = '):
                    if jira == None:
                        print "Using default value from template file", line
                        outfile.write(line)
                        continue
                    else:
                        line = line[0:len('JIRA_ID = ')] + jira + '\n'
                        outfile.write(line)
                elif line.startswith('EVENT_TAG = '):
                    if event == None:
                        print "Using default value from template file", line
                        outfile.write(line)
                        continue
                    else:
                        line = line[0:len('EVENT_TAG = ')] + event + '\n'
                        outfile.write(line)
                elif line.startswith('RUN_COMMENT = '):
                    if run_comment == None:
                        print "Using default value from template file", line
                        outfile.write(line)
                        continue
                    else:
                        line = line[0:len('RUN_COMMENT = ')] + run_comment + '\n'
                        outfile.write(line)
                elif line.startswith('NITE = '):
                    if night == None:
                       print "Must provide a night to process"
                       sys.exit(0)
                    else: 
                        line = line[0:len('NITE = ')] + night + '\n'
                        outfile.write(line)
                elif line.startswith('PROJECT = '):
                    if project == None:
                        print "Using ACT project..."
                        line = line[0:len('PROJECT = ')] + 'ACT\n'
                        outfile.write(line)
                    else:
                        line = line[0:len('PROJECT = ')] + project + '\n'
                        outfile.write(line)
                elif line.startswith('biascor_run = '):
                    if biascor == None:
                        print "Using default value from template file", line
                        outfile.write(line)
                        continue
                    else:
                        line = line[0:len('biascor_run = ')] + biascor + '\n'
                        outfile.write(line)
                elif line.startswith('flatcor_run = '):
                    if flatcor == None:
                        print "Using default value from template file", line
                        outfile.write(line)
                        continue
                    else:
                        line = line[0:len('flatcor_run = ')] + flatcor + '\n'
                        outfile.write(line)
                else:
                    outfile.write(line)
                    
                    
    if os.path.isfile(submit_file):
        #if the file exists, check if it is empty
        try:
            s = os.stat(submit_file)
            if s.st_size == 0:
                print "The submit file {} is empty".format(submit_file)
                return sys.exit(1)
        except OSError as e:
            print e
            return sys.exit(2)
        
        #File exists and is not empty, so return 1 (True)
        return submit_file
    else:
        #file does not exists, return 0
        print "File doesn't exists"
        sys.exit(1)
        
    
    


def submit_diffim(red_run, field, jira_id):
    """
    Reads submit file for difference image from submit directory and updates 
    parameters
    :red_run : red run from firstcut_snese
    :field : field to do diffim that is present in red_run 
    :param jira_id: Optional Jira_id where notes are taken
        
    """
    
    
def desstat():
    """
    Check runs running and look for how many nodes are available to submit a new run
    """

    try:
        des_home = os.environ['DES_HOME']
        print "DES_HOME is define. Continue processig"
    except KeyError:
        des_home  = None
    
    if des_home is None:
        print "setting up desdmsoft current"
        set_deshome = subprocess.Popen('$EUPS_DIR/bin/eups_setup desdmsoft 7.2.10+1', shell=True)
        out, err = set_deshome.communicate()
        print "Return code: ", set_deshome.returncode
        #print out.rstrip(), err.rstrip()

        #os.system('/home/ricardoc/DESDM/Y1-SNe/scripts/setup_desdmsoft.csh')

    
    snese_path = '/home/ricardoc/DESDM/Y1-SNe/snese/'
    run_running = os.path.join(snese_path,'temp_desstat.txt')
    
    print run_running
    
    p = subprocess.Popen('desstat', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    
    with open(run_running,'w') as outfile:
        #subprocess.call(['desstat'], stdout = outfile)
        for line in p.stdout.readlines():
            outfile.write(line)

    #check the temp desstat file
    #ID T PRJ PIPELINE RUN                         BLOCK                          SUBBLOCK        STAT              OPERATOR RUNSITE
    #============================================================================================================================================
    #35744   OPS snese    20130912072031_20130911     make_craymask                  pipelines       (0/0/2/2)         ricardoc iforge
    #35848   OPS firstcut 20130909091938_20130908     endrun                         runpost         RUN               mgelman  desjobc3
    #35853   OPS firstcut 20130911103044_20130910     endrun                         runpost         RUN               mgelman  desjobc3

    #check if the file is empty
    try:
        s = os.stat(run_running)
        if s.st_size == 0:
            print "The desstat file {} is empty".format(run_running)
    except OSError as e:
        print e
        sys.exit(1)

    #file is not empty. Then read it.
    id = []
    project = []
    pipeline = []
    run = []
    block = []
    subblock = []
    stat = []
    operator= []
    runsite = []
    
    if os.path.isfile(run_running):
        #with open('/home/ricardoc/DESDM/Y1-SNe/snese/temp','r') as infile:
        with open(run_running,'r') as infile:
            for line in infile:
                #print " line is ", line.rstrip() 
                if line in ['\n', '\r\n']:
                    continue
                elif line.startswith('No DES'):
                    print "No DES Jobs running. Code is able to submit a job"
                    stat = ['NO_DES_JOBS']
                    continue
                elif ('ID' in line) or (line.startswith('=')):
                        #print line
                        continue
                else:
                    line = line.rstrip()
                    data = line.split()
                    #print line
                    #print data
                    if len(data) == 9:
                        id.append(data[0])
                        project.append(data[1])
                        pipeline.append(data[2])
                        run.append(data[3])
                        block.append(data[4])
                        subblock.append(data[5])
                        stat.append(data[6])
                        operator.append(data[7])
                        runsite.append(data[8])
                    elif len(data) == 8:
                        continue
                    else:
                        print "Missing column in destest output file. Can't determine exactly how many processes are running"
                        print "check file {} ".format(run_running)
                        sys.exit(1)

    
    
    #print RUNSITE,STAT
    
    return stat
                    

    
def check_total_nodes_available(status_runs):
    """
    Check The total number of cores available on iforge.
    If the total number of running of core currently used is 9, then code
    can submit the current run.
    If not, then it has to wait until condition is true
    
    :param status_runs: List with all processing running (output from desstat())
    """
        
    #initialize current "queue, running and halt" proccessing
    total_halt = 0
    total_used_nodes = 0
    total_queue = 0
    
    #check using pipeline, stat and block
    for value in status_runs: 
        if value.startswith('('):
            halt, queue, running, total = value.split('/')
            halt = halt.lstrip('(')
            total = total.rstrip(')')
            total_used_nodes = total_used_nodes + int(running)
            total_queue = total_queue + int(queue)
            total_halt = total_halt + int(halt)

    print "Total Number of cores currently being used: ", total_used_nodes
    
    return total_used_nodes
    
def submit_run(dessubmit_file):
    """
    Submit run for processing.
    
    :param dessubmit_file: submit file created automatically for current run to process
    
    """
    
    try:
        des_home = os.environ['DES_HOME']
        print "DES_HOME is define. Continue processig"
    except KeyError:
        des_home  = None

    if des_home is None:
        print "setting up desdmsoft current"
        set_deshome = subprocess.Popen('$EUPS_DIR/bin/eups_setup desdmsoft 7.2.10+1', shell=True)
        #out, err = set_deshome.communicate()
        print "Return code: ", set_deshome.returncode
        #print out.rstrip(), err.rstrip()

        #os.system('/home/ricardoc/DESDM/Y1-SNe/scripts/setup_desdmsoft.csh')

    
    #NOTES:
    #I NEED TO CHANGE TO THE SUBMIT DIRECTORY BEFORE I TRY TO SUBMIT FILE
        
    dessubmit_path_log = '/home/ricardoc/DESDM/Y1-SNe/log/'
    
    #extract the filename from the dessubmit_file. 
    submit_path, submit_file = os.path.split(dessubmit_file)
    submit_file_log = submit_file.replace('des','log')
    
    full_path_submit_file_log = os.path.join(dessubmit_path_log,submit_file_log)
    
    #move to submit directory
    os.chdir(submit_path)
    
    print "submit_file log is", full_path_submit_file_log
    
    dessubmit_cmd = 'dessubmit ' + submit_file
    
    print "dessubmit cmd", dessubmit_cmd
    #submit the process and create the log file
    p = subprocess.Popen(dessubmit_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print "Return code: ", p.returncode
    
    with open(full_path_submit_file_log,'w') as outfile:
        for line in p.stdout.readlines():
            outfile.write(line)




def main():

    
    #list all observation type present in OBS_SET_TYPE table.
    if list_type:
        mycur = connectDB("db-desoper")
        list_observations_types(mycur)
  

    #make Database connection
    cursor_connect = connectDB("db-desoper")
    
    #check if the exposure in obs_set_exposure (written in from manifest) are already
    #ingested in DB.
    queryitems = ["ose.exposurename", "os.type", "e.exposurename"]
    #queryitems = ["exposurename, nite, band, field, sn_sequence_status"]
    querylist = ",".join(queryitems)

    #If no nite is given as input, then assume that we are processing data
    #from current observed night.
    
    #Begin assuming the exposure in manifest file have not been ingested to DB
    are_all_ingested = False
    
    while not are_all_ingested:
        exposure_manifest, ingested_exposures, type = check_exposures_in_manifest(night,querylist,cursor_connect)
        #check if all exposures in manifest have been ingested.
        are_all_ingested = manifest_exposures_equal_to_ingested_exposures(exposure_manifest,ingested_exposures,verbose=verbose)
        
        if are_all_ingested:
            print "All Exposures have arrived and ingested to NCSA"
        else:
            print "Waiting 60 sec for cheking again DB if files has arrived"
            time.sleep(60)

    if len(type) > 0:
        exposure_type = type[0]
    print "creating idlist for exposures to process"        
    success_create_list, idlist_filename = create_idlist_for_submit_file(ingested_exposures,night,cursor_connect,exposure_type,verbose=verbose)

    #create snese submit file
    if success_create_list:
        print "creating snese submit file"
        if night == None:
            print "Need a night to process"
        
        #if event == None:
        exposure_type_tag = exposure_type.replace('-','_')
        event = 'snese_' + night + '_' + exposure_type_tag + '_attempt1'
        #if run_comment == None:
        run_comment = 'snese_' + night + '_' + exposure_type_tag + '_attempt1'
        
        print "exp_type is", exposure_type
        #create submit file. If creation is successful, then return True
        #if file created is empty or failed, then return False
        created_submit_file = create_submit_snese(night=night, jira=None, event=event, project=project, run_comment = run_comment, idfile = idlist_filename, type = exposure_type, biascor = biascor, flatcor = flatcor)

        print "created submit file", created_submit_file
        
        #before submit run, need to check if there are nodes available.        
        stat = desstat()
        
        #check of many nodes are currently being used on iforge
        total_used_nodes = check_total_nodes_available(stat)
        
        #cmd = "ls -l ./"
        #p = subprocess.Popen(cmd , shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #out, err = p.communicate()
        #print "Return code: ", p.returncode
        #print "xxxx",out.rstrip(), err.rstrip()
        
        #submit run only if the total number of running run is less than 9 
        #dessubmit = 'dessubmit ' + created_submit_file
        
        submit_run(created_submit_file)
        
        #if (total_running < 9):
            #p.subprocess.Popen('dessubmit %s' % (created_submit_file), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                
        

if  __name__ == '__main__':
    main()


