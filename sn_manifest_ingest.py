#!/usr/bin/env python
"""Select a set of raw observations from the DB and populates 
a test table (SN_MANIFEST)

.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>



Given a field and a night, it finds all the src images in the DB and populates
the SN_MANIFEST table

Some notes from Ankit.
I need to create two tables. One that contains all the observations for a given night
and anothe that contains the status of observations for a given SN field.

Second table would be:



MANIFEST
MANIFEST_ID             NOT NULL NUMBER(22)  
SEQUENCE_ID            NOT NULL NUMBER(22)
EXPOSURENAME          VARCHAR2(100)
INSERT_TIME               VARCHAR(80)
 

SEQUENCE_OBSERVATIONS
SEQUENCE_ID                          NOT NULL NUMBER(22)
NUMBER_OBSERVATIONS       NUMBER(22)
SEQUENCE_TYPE                     VARCHAR(80)
NITE                                   VARCHAR(80)

pre-defined sequence types:

SN-X3-g
SN-X3-r
SN-X3-i
SN-X3-z
SN-C3-g
SN-C3-r
SN-C3-i
SN-C3-z
SN-C1
SN-C2
SN-S1
SN-S2
SN-E1
SN-E2
SN-X1
SN-X2
zero
domeflat-u
domeflat-g
domeflat-r
domeflat-i
domeflat-z
domeflat-Y
focus
junk
test
survey
standards
unknown

This table is in ricardoc user account..

usage: sn_manifest_ingest.py [-h] -n NIGHT [-f FIELD] [-v]

Given a night and optionally a field name for SNe, it populates sn_manifest
table

optional arguments:
  -h, --help            show this help message and exit
  -n NIGHT, --nite NIGHT
                        Single Night
  -f FIELD, --field FIELD
                        Optional field name
  -v, --verbose         Use --verbose for verbose output

"""

# imports
import coreutils.desdbi
import argparse
import os
import time
import string

parser = argparse.ArgumentParser(description="Given a night and optionally a field name for SNe, it populates sn_manifest table")
parser.add_argument("-n", "--nite", required=True, dest="night", help="Single Night")
parser.add_argument("-f", "--field", required=False, dest="field", help="Optional field name")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")

args = parser.parse_args()
night = args.night
field = args.field
verbose = args.verbose


#    Turn on svn keyword substitution for Revision (above) with:
#    bash$ svn propset svn:keywords "Revision" foo_quality.py
__version__ = "$Revision: 187 $"

def main():
    
    try:
        desdmfile = os.environ["des_services"]
    except KeyError:
        desdmfile = None
    dbh = coreutils.desdbi.DesDbi(desdmfile,"db-desoper")
    cur = dbh.cursor()



    queryitems = ["exposurename, nite, band, object"]
    querylist = ",".join(queryitems)

        
    #if verbose:
    #    print "Using a search box of dRA %.2f degrees, dDec %.2f degrees\n" % (dra_in_deg,ddec_in_deg)

    #I used this to open a file and read input from it
    #allRows =[]
    #with open(file,'r') as f:
    #    for line in f:
    #        line = line.rstrip()
    #        allRows.append(line)

    if field:
        my_field = 'SN-' + field
        print my_field

    if field is None:
        query = """select distinct %s
        from exposure
        where nite = '%s' 
        and object like '%%SN-%%' 
        and exposuretype='src' order by exposurename""" % (querylist,  night )
        print query
    else:
        query = """select distinct %s
        from exposure
        where nite = '%s'
        and object like '%%%s%%' 
        and exposuretype='src' order by exposurename""" % ( querylist, night, my_field )

        
    cur.arraysize = 1000 # get 1000 at a time when fetching
    if verbose:
        print "Performing query: %s \n\n" % query
    cur.execute(query)
    row = cur.fetchall()
    

    #counter to simulate when the observations are done
    count=0
    for line in row:
        #print row[count]
        exposurename = line[0]
        nite = line[1]
        band = line[2]
        object = line[3]
        curr_object = object
        curr_field_pos = string.find(curr_object,'SN-')
        curr_field = curr_object[curr_field_pos:curr_field_pos+5]
        if count==0:
            print 'current field observing is %s %s' % (curr_field,line[0])
        
        #I need to wait 5 minutes before I submit the next processing.
        #cur.execute(sqlInsert, line)
        #cur.commit()
        
        if count == len(row)-1:
            next_field = curr_field
        else:
            next_row = row[count+1]
            #print prev_row
            next_object = next_row[3]
            next_field_pos = string.find(next_object,'SN-')
            next_field = next_object[next_field_pos:next_field_pos+5]
            
        #print 'previous field is %s and curr_file %s' % (next_field, curr_field)


        insertArray = []
        #now check if field observations is finished and add the corresponding trigger to the table
        if (curr_field != next_field) or (count == len(row)-1):
            sn_seq_status = 1  #Observations for field have not finished
            insertArray.append([exposurename, nite, band, curr_field, object, sn_seq_status])
            sqlInsert = 'INSERT INTO SN_MANIFEST VALUES (sn_manifest_seq.nextval, :1, :2, :3, :4, :5, :6)'
        else:
            sn_seq_status =0 #
            insertArray.append([exposurename, nite, band, curr_field, object, sn_seq_status])
            sqlInsert = 'INSERT INTO SN_MANIFEST VALUES (sn_manifest_seq.nextval, :1, :2, :3, :4, :5, :6)' 
        #print sqlInsert
        
 
        #now let's populate the table every 90 sec
        #Code is sleeping for 90 secs
        if count > 0:
            time.sleep(90)
        
        cur.executemany(sqlInsert,insertArray)
        dbh.commit()
 
        count += 1
        
        print "Inserted %s of band: %s , nite:  %s , object name:  %s \n" % (exposurename, band, nite, object)

if  __name__ == '__main__':
    main()

