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
         print "error: could not change to "+path
         sys.exit("ending session")

    ftp.retrbinary("RETR "+file_name, open(os.path.join(destination,file_name),"wb").write)
    print file_name + " downloaded"



def download_whole_directory(ftp, path, destination, recursive=False):
    print path
    try:
        ftp.cwd(path)
        os.chdir(destination)
    except OSError:
        pass
    except ftplib.error_perm:
        print "error: could not change to "+path
        sys.exit("ending session") 

    filelist=ftp.nlst()
    print "files:"
    print filelist
    for file in filelist:
        print "downloading...", os.path.join(path,file)
        try:
            ftp.cwd(os.path.join(path,file)+"/")

            if recursive == True:
                os.mkdir(os.path.join(destination,file))              
                download_whole_directory(ftp, path+file+"/", os.path.join(destination,file))
        except ftplib.error_perm:
            ftp.cwd(path)
            print "downloading", file
            #print ftp.nlst()
            os.chdir(destination)
            try:
                ftp.retrbinary("RETR "+file, open(file, "wb").write)
                print file + " downloaded"
            except ftplib.error_perm:
                print 'could not download file/dir: %s' % file