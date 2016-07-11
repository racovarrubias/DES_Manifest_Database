#! /usr/bin/env python


def split_section(section):
#
#   Function to split a section keyword string '[x1:x2,y1:y2]' into its component values.
#
	xy_range=section[1:-1].split(',')
        x_range=xy_range[0].split(':')
        y_range=xy_range[1].split(':')
	section_val={}
	section_val["x1"]=int(x_range[0])
	section_val["x2"]=int(x_range[1])
	section_val["y1"]=int(y_range[0])
	section_val["y2"]=int(y_range[1])
	return section_val
	

##########################################################################################
def quickimstat2(indata,section,clipsig,maxiter,converge_num,verbose):
#
#   Function that reforms a subsection of science data array and obtains statistics
#
        tmp_data=indata[section["y1"]:section["y2"],section["x1"]:section["x2"]]
        tmp_data2=tmp_data.reshape(tmp_data.size,)
        retval=medclip(tmp_data2,clipsig,maxiter,converge_num,verbose)
        avgval=retval[0]; medval=retval[1]; stdval=retval[2]
#       prf = 'MEANCLIP:'
#       print '%s %.1f-sigma clipped mean' % (prf, clipsig)
#       print '%s Mean computed in %i iterations' % (prf, iter)
#       print '%s Mean = %.6f, sigma = %.6f' % (prf, avgval, stdval)
        return avgval,medval,stdval

##########################################################################################
def quickimstat_jump(indata,sectiona,sectionb,clipsig,maxiter,converge_num,verbose):
	tmp_data_jump=numpy.zeros((sectiona["y2"]-sectiona["y1"]),dtype=numpy.float32)
	for iy in range(sectiona["y1"],sectiona["y2"]):
		maxiter=3
		converge_num=0.1
		clipsig=3.0
	        tmp_data_a=indata[iy:iy+1,sectiona["x1"]:sectiona["x2"]]
	        tmp_data_a2=tmp_data_a.reshape(tmp_data_a.size,)
	        tmp_data_b=indata[iy:iy+1,sectionb["x1"]:sectionb["x2"]]
	        tmp_data_b2=tmp_data_b.reshape(tmp_data_b.size,)
	        retvala=medclip(tmp_data_a2,clipsig,maxiter,converge_num,verbose)
	        retvalb=medclip(tmp_data_b2,clipsig,maxiter,converge_num,verbose)
		temp_comp=retvala[0]/retvalb[0]
		if (temp_comp < 0.):
			tmp_data_jump[(iy-sectiona["y1"])]=0.0
		else:
			tmp_data_jump[(iy-sectiona["y1"])]=retvala[0]/retvalb[0]
#		if (iy%100 == 0):
#			print "###############################"
#			print sectiona
#			print sectionb
#			print iy,indata[iy,1023],indata[iy,1024]
#			print iy,indata[iy:iy+1,1023:1024],indata[iy:iy+1,1024:1025]
#			print "  a2: ",tmp_data_a2
#			print "  b2: ",tmp_data_b2
#			print " ret: ",iy,retvala,retvalb

	tmp_data_jump2=numpy.sqrt(tmp_data_jump)
	
	maxiter=10
	converge_num=0.0001
	clipsig=3.0
        retval_jump=medclip(tmp_data_jump2,clipsig,maxiter,converge_num,verbose)
        avgval=retval_jump[0]; medval=retval_jump[1]; stdval=retval_jump[2]
#       prf = 'MEANCLIP:'
#       print '%s %.1f-sigma clipped mean' % (prf, clipsig)
#       print '%s Mean computed in %i iterations' % (prf, iter)
#       print '%s Mean = %.6f, sigma = %.6f' % (prf, avgval, stdval)
        return avgval,medval,stdval

##########################################################################################
def medclip(data,clipsig,maxiter,converge_num,verbose):
#
#   Function to provide statistics over a vector
#
        ct = data.size
        iter = 0; c1 = 1.0 ; c2 = 0.0
        while (c1 >= c2) and (iter < maxiter):
                lastct = ct
                medval = numpy.median(data)
                sig = numpy.std(data)
                wsm = numpy.where( abs(data-medval) < clipsig*sig )
                ct = len(wsm[0])
                if ct > 0:
                        data = data[wsm]

                c1 = abs(ct - lastct)
                c2 = converge_num * lastct
#               print iter,len(data)
                iter += 1
#       End of while loop

        avgval = numpy.mean(data)
        stdval = numpy.std(data)
        return avgval,medval,stdval

if __name__ == "__main__":

    # Importing necessary modules
    from optparse import OptionParser
    import os
    import coreutils.desdbi
    import fnmatch
    import tempfile
    import re
    import time
    import sys
    import numpy
    import fitsio
#    import pyfits
#    from numpy import * 

    svnid="$Id: gather_flatframe_qa.py 16000 2014-01-06 21:00:00Z rgruendl $"

#   Creating command-line options
    parser = OptionParser(usage=__doc__)
    parser.add_option("-e","--explist", dest="explist", default=None, help="Required (comma separated list of exposure numbers)")
    parser.add_option("-p","--path", dest="path", default=".", help="Required path where temporary crosstalk products are waiting (default=.)")
    parser.add_option("-c","--ccd", dest="ccd_restrict", default="All", help="Restrict to specific CCDs (default=All)")
    parser.add_option("-u","--updateDB", action="store_true", dest="updateDB", default=False,
			help="Flag for program to DIRECTLY update DB (firstcut_eval).")
    parser.add_option("-D","--DB_file", action="store_true", dest="DB_file", default=False,
			help="Flag for optional output of DB update file")
    parser.add_option("-O","--Override", action="store_true", dest="override", default=False, help="Flag to override non-standard files")
    parser.add_option("-f","--froot", dest="froot", default="temp", help="Root for output file names (default=temp)")
    parser.add_option("-v","--verbose",action="store_true",dest="verbose",default=False,
			help="Print queries and show progress")
    parser.add_option("-s","--section", dest="section", default="db-desoper", 
			help="Section of .desservices file with connection into (default is db-desoper)")
    (opts,args) = parser.parse_args()

    if (opts.explist is None):
	parser.print_help()
	exit(1)	
    else:
	prelim_expnum_list=[]
	tmp_list=opts.explist.split(',')
	for tmp_entry in tmp_list:
		if (tmp_entry.strip() != ''):
			prelim_expnum_list.append(int(tmp_entry.strip()))
	for expnum_index, expnum_val in enumerate(prelim_expnum_list):
		if (expnum_index == 0):
			exp_constraint=' and e.expnum in ('+'%d' % expnum_val
		else:
			exp_constraint=exp_constraint+','+'%d' % expnum_val
	exp_constraint=exp_constraint+')'
	print exp_constraint

    verbose = opts.verbose
    ccd_list=[]
    if (opts.ccd_restrict == "All"):
	for ccd in range(1,63):
		ccd_list.append(int(ccd))
    else:
	ccd_list=opts.ccd_restrict.split(',')
    print ccd_list
#
#   Sigma clipping algorithm parameters
#
    maxiter=10
    converge_num=0.0001
    clipsig=3.0
    jumpwid=25

    try:
	desdmfile = os.environ["des_services"]
    except KeyError:
	desdmfile = None
    dbh = coreutils.desdbi.DesDbi(desdmfile,opts.section)
    cur = dbh.cursor()

    old_db_table="gruendl.flat_frame_qa"
    new_db_table="gruendl.flat_frame_qa_v2"

#
#   First query the database for information usually found when working on "raw" products.
#
    queryitems = ["e.id", "e.expnum", "e.band", "e.exptime", "e.object", "e.obstype" ]
    coldict={}
    for index, item in enumerate(queryitems):
        coldict[item]=index
    querylist = ",".join(queryitems)
    query = """select %s from exposure e where e.project='DTS' %s """ % (querylist,exp_constraint)

    if opts.verbose:
	print query
    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)

    file_list=[]
    for item in cur:
	tmp_dict={}
	tmp_dict["expnum"]=int(item[coldict["e.expnum"]])
	tmp_dict["eid"]=int(item[coldict["e.id"]])
	tmp_dict["band"]=item[coldict["e.band"]]
	tmp_dict["exptime"]=float(item[coldict["e.exptime"]])
	tmp_dict["obstype"]=item[coldict["e.obstype"]]
	tmp_dict["object"]=item[coldict["e.object"]]
	if ((tmp_dict["band"] in ["u","g","r","i","z","Y"])and(tmp_dict["obstype"]=="dome flat")):
		file_list.append(tmp_dict)
		print("  Process: {:8d} {:1s} {:12s} {:6.1f} {:s} ".format(tmp_dict["expnum"],tmp_dict["band"],tmp_dict["obstype"],tmp_dict["exptime"],tmp_dict["object"]))
	elif (opts.override):
		file_list.append(tmp_dict)
		print(" Override: {:8d} {:1s} {:12s} {:6.1f} {:s} ".format(tmp_dict["expnum"],tmp_dict["band"],tmp_dict["obstype"],tmp_dict["exptime"],tmp_dict["object"]))
	else:
		print(" Skipping: {:8d} {:1s} {:12s} {:6.1f} {:s} ".format(tmp_dict["expnum"],tmp_dict["band"],tmp_dict["obstype"],tmp_dict["exptime"],tmp_dict["object"]))

    print len(file_list)

#    exit(0)

    if (opts.DB_file):
	fdbout=open("%s.dbupdate"% (opts.froot), 'w')

    run='fromsrc'
    icount=0
    t0=time.time()
    for file_record in file_list:
	for ccd in ccd_list:
		fname=os.path.join(opts.path,("DECam_%08d_%02d.fits" % (file_record["expnum"],ccd)))
		if (not(os.path.isfile(fname))):
			print "# Warning: ",fname," was not found."
		else:
			icount=icount+1
#
#			Now check to make see whether an entry already exists in the OLD or NEW database table 
#			For the OLD table query go ahead and get the values so there is no need to recompute.
#
			queryitems = ["f.expnum", "f.ccd", "f.band", "f.source", "f.med_a", "f.med_ar", "f.med_aw", "f.med_b", "f.med_br", "f.med_bw", "f.rms_a", "f.rms_ar", "f.rms_aw", "f.rms_b", "f.rms_br", "f.rms_bw"] 
			coldict={}
			for index, item in enumerate(queryitems):
				coldict[item]=index
			querylist = ",".join(queryitems)
			query = """select %s from %s f where f.source='%s' and f.expnum=%d and f.ccd=%d """ % (querylist,old_db_table,run,file_record["expnum"],ccd)
			if opts.verbose:
				print query
			cur.arraysize = 1000 # get 1000 at a time when fetching
			cur.execute(query)

			old_DBentry=False
			for item in cur:
				if (item[0]>0):
					old_DBentry=True
					statval_a=[]
					statval_a.append(0.0)
					statval_a.append(float(item[coldict["f.med_a"]]))
					statval_a.append(float(item[coldict["f.rms_a"]]))
					statval_ar=[]
					statval_ar.append(0.0)
					statval_ar.append(float(item[coldict["f.med_ar"]]))
					statval_ar.append(float(item[coldict["f.rms_ar"]]))
					statval_aw=[]
					statval_aw.append(0.0)
					statval_aw.append(float(item[coldict["f.med_aw"]]))
					statval_aw.append(float(item[coldict["f.rms_aw"]]))
					statval_b=[]
					statval_b.append(0.0)
					statval_b.append(float(item[coldict["f.med_b"]]))
					statval_b.append(float(item[coldict["f.rms_b"]]))
					statval_br=[]
					statval_br.append(0.0)
					statval_br.append(float(item[coldict["f.med_br"]]))
					statval_br.append(float(item[coldict["f.rms_br"]]))
					statval_bw=[]
					statval_bw.append(0.0)
					statval_bw.append(float(item[coldict["f.med_bw"]]))
					statval_bw.append(float(item[coldict["f.rms_bw"]]))


			query = """select count(*) from %s f where f.source='%s' and f.expnum=%d and f.ccd=%d """ % (new_db_table,run,file_record["expnum"],ccd)
			if opts.verbose:
				print query
			cur.arraysize = 1000 # get 1000 at a time when fetching
			cur.execute(query)

			new_DBentry=False
			for item in cur:
				if (item[0]>0):
					new_DBentry=True

#
#			If a NEW database entry already exists then skip this exposure/ccd
#
			if (new_DBentry):
				print "# Note: Entry present in %s for expnum=%d, ccd=%02d " % (new_db_table,file_record["expnum"],ccd)
			else:		
#
#				If file is (un)compressed then need to look at different HDUs
#
				if (fname[-2:] == "fz"):
					sci_hdu, msk_hdu, wgt_hdu = (1,2,3) # for .fz
				else:
					sci_hdu, msk_hdu, wgt_hdu = (0,1,2) # for .fits (or .gz)

				fits = fitsio.FITS(fname,'r') # Change to 'rw' if you want
				h = fits[sci_hdu].read_header()
				SCI=fits[sci_hdu].read()
				corn_ampa=split_section(h["AMPSECA"])
				ccdnum=int(h["CCDNUM"])

				if (corn_ampa["x1"]==1):
					asection={"x1":256,"x2":768,"y1":1024,"y2":3072}
					asection_r={"x1":256,"x2":768,"y1":50,"y2":100}
					asection_w={"x1":256,"x2":768,"y1":3996,"y2":4046}
					bsection={"x1":1280,"x2":1792,"y1":1024,"y2":3072}
					bsection_r={"x1":1280,"x2":1792,"y1":50,"y2":100}
					bsection_w={"x1":1280,"x2":1792,"y1":3996,"y2":4046}
					asection_j={"x1":int(1024-jumpwid),"x2":1024,"y1":50,"y2":4046}
					bsection_j={"x1":1024,"x2":int(1024+jumpwid),"y1":50,"y2":4046}
				else:
					asection={"x1":1280,"x2":1792,"y1":1024,"y2":3072}
					asection_r={"x1":1280,"x2":1792,"y1":3996,"y2":4046}
					asection_w={"x1":1280,"x2":1792,"y1":50,"y2":100}
					bsection={"x1":256,"x2":768,"y1":1024,"y2":3072}
					bsection_r={"x1":256,"x2":768,"y1":3996,"y2":4046}
					bsection_w={"x1":256,"x2":768,"y1":50,"y2":100}
					asection_j={"x1":1024,"x2":int(1024+jumpwid),"y1":50,"y2":4046}
					bsection_j={"x1":int(1024-jumpwid),"x2":1024,"y1":50,"y2":4046}

#	print(" {:2d}  {:4d} {:4d} {:4d} {:4d}  {:4d} ".format(ccd,corn_ampa["x1"],corn_ampa["x2"],corn_ampa["y1"],corn_ampa["y2"],asection_r["y1"]))
#	print ccdnum,ccd,corn_ampa["x1"],corn_ampa["y1"],asection,asection_r,bsection,bsection_r

				if (not(old_DBentry)):
					maxiter=10
					converge_num=0.0001
					clipsig=3.0
	
					statval_a=quickimstat2(SCI,asection,clipsig,maxiter,converge_num,opts.verbose)
					statval_ar=quickimstat2(SCI,asection_r,clipsig,maxiter,converge_num,opts.verbose)
					statval_aw=quickimstat2(SCI,asection_w,clipsig,maxiter,converge_num,opts.verbose)

					statval_b=quickimstat2(SCI,bsection,clipsig,maxiter,converge_num,opts.verbose)
					statval_br=quickimstat2(SCI,bsection_r,clipsig,maxiter,converge_num,opts.verbose)
					statval_bw=quickimstat2(SCI,bsection_w,clipsig,maxiter,converge_num,opts.verbose)

				maxiter=2
				converge_num=0.1
				clipsig=3.0
				statval_abjump=quickimstat_jump(SCI,asection_j,bsection_j,clipsig,maxiter,converge_num,opts.verbose)

				fits.close()

				if ((opts.updateDB)or(opts.DB_file)):
					if (not(old_DBentry)):
						db_cname="INSERT INTO %s (EXPNUM,CCD,BAND,SOURCE,MED_A,MED_AR,MED_AW,MED_B,MED_BR,MED_BW,RMS_A,RMS_AR,RMS_AW,RMS_B,RMS_BR,RMS_BW) " % ( old_db_table )
					        db_value=""" VALUES (%d, %d, '%s', '%s', %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f ) """ % (file_record["expnum"],ccd,file_record["band"],run,statval_a[1],statval_ar[1],statval_aw[1],statval_b[1],statval_br[1],statval_bw[1],statval_a[2],statval_ar[2],statval_aw[2],statval_b[2],statval_br[2],statval_bw[2])
						insert_command=db_cname+db_value
						if (opts.updateDB):
							cur.execute(insert_command)
						if (opts.DB_file):
							fdbout.write("{:s};\n".format(insert_command))
						else:
							print file_record["expnum"],ccd,file_record["band"],statval_a[1],statval_a[2],statval_b[1],statval_b[2],statval_abjump[1],statval_abjump[2]


					if (not(new_DBentry)):
						db_cname="INSERT INTO %s (EXPNUM,CCD,BAND,SOURCE,MED_A,MED_AR,MED_AW,MED_B,MED_BR,MED_BW,MED_AB,RMS_A,RMS_AR,RMS_AW,RMS_B,RMS_BR,RMS_BW,RMS_AB) " % ( new_db_table )
					        db_value=""" VALUES (%d, %d, '%s', '%s', %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.5f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.5f) """ % (file_record["expnum"],ccd,file_record["band"],run,statval_a[1],statval_ar[1],statval_aw[1],statval_b[1],statval_br[1],statval_bw[1],statval_abjump[1],statval_a[2],statval_ar[2],statval_aw[2],statval_b[2],statval_br[2],statval_bw[2],statval_abjump[2])
						insert_command=db_cname+db_value

						if (opts.updateDB):
							cur.execute(insert_command)
						if (opts.DB_file):
							fdbout.write("{:s};\n".format(insert_command))
						else:
							print file_record["expnum"],ccd,file_record["band"],statval_a[1],statval_a[2],statval_b[1],statval_b[2],statval_abjump[1],statval_abjump[2]
			if (icount%100 == 0):
				t1=time.time()	
				print("@ Completed {:6d} images in {:8.1f} seconds ({:8.2f} seconds/image).".format(icount,(t1-t0),((t1-t0)/icount)))
	
    if(opts.updateDB):
	dbh.commit()
	print "DB update complete and committed"
    if (opts.DB_file):
	fdbout.close()


