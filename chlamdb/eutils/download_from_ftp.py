#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

# source ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

#import tarfile
#import re
#from ftplib import FTP
#import operator
import sys
import ftplib
import os




def download_one_file(ftp, path, destination, file_name):
    try: 
        ftp.cwd(path)
        os.chdir(destination)
    except OSError:
        pass
    except ftplib.error_perm:
         print ("error: could not change to "+path)
         sys.exit("ending session")

    ftp.retrbinary("RETR "+file_name, open(os.path.join(destination,file_name),"wb").write)
    print (file_name + " downloaded")



def download_whole_directory(ftp, path,
                             destination,
                             recursive=False,
                             only_gbk=False,
                             only_fna=False):
    #print path
    #print 'recursive:', recursive
    try:
        ftp.cwd(path)
        os.chdir(destination)
    except OSError:
        print ('could not reach directory: %s' % path)
    except ftplib.error_perm:
        print ("error: could not change to "+path)
        sys.exit("ending session") 

    filelist=ftp.nlst()

    if only_fna and only_gbk:
        raise('choose either gbk or fna')

    if only_gbk:
        filelist = [i for i in filelist if 'genomic.gbff.gz' in i]
    if only_fna:
        filelist = [i for i in filelist if 'genomic.fna.gz' in i]

    #print "files:"
    #print filelist

    if filelist[0] == 'assembly_status.txt':
        return False
    
    for file in filelist:
        #print 'dir:', ftp.pwd()
        #print "downloading...", os.path.join(path,file)
        if recursive == True:
            try:
                ftp.cwd(os.path.join(path, file)+"/")
                os.mkdir(os.path.join(destination,file))
                download_whole_directory(ftp, path+file+"/", os.path.join(destination, file))

            except ftplib.error_perm:
                #print "downloading", file
                os.chdir( destination)
                try:
                    ftp.retrbinary("RETR "+file, open(file, "wb").write)
                    #print file + " downloaded"
                except ftplib.error_perm:
                    print (ftp.nlst())
                    print (ftp.pwd())
                    print ('could not download file/dir: %s' % file)
        else:
                #print "downloading", file
                os.chdir(destination)
                try:
                    ftp.retrbinary("RETR "+file, open(file, "wb").write)
                    #print file + " downloaded"
                except ftplib.error_perm:
                    print (ftp.nlst())
                    print (ftp.pwd())
                    print ('could not download file/dir: %s' % file)
