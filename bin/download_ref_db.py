#!/usr/bin/env python


import ftplib
import re
import io
import os
import threading
import time
import argparse
import gzip


def download_refseq(download_dir, n_retry=5):
    ftp = ftplib.FTP("ftp.ncbi.nih.gov")
    ftp.login("anonymous")
    ftp.cwd("/refseq/release/complete")

    # check if there are any pre-existing files in the download directory
    # to avoid downloading them a second time
    # NOTE: this assumes that the files are either completely downloaded
    # or not present (i.e. removed from the disk in case of problem)
    existing_files = set()
    for current_file in os.listdir(download_dir):
        if os.path.isfile(download_dir+"/"+current_file):
            existing_files.add(current_file)

    # it would probably be better to use startsWith and check the last
    # letters to make sure it ends with faa.gz, instead of using regular 
    # expression
    nr_re = re.compile("complete.nonredundant_protein.(\d)*.protein.faa.gz")
    nr_filelist = (i for i in ftp.nlst() if not re.match(nr_re, i) is None and
            not i in existing_files)

    for f in nr_filelist:
        failed = 0
        complete = False
        while not complete:
            output_file = open(download_dir+"/"+f, "wb")
            try:
                ftp.retrbinary("RETR "+f, output_file.write)
                output_file.close()
                complete = True
            except:
                output_file.close()
                failed += 1
                if failed==n_retry:
                    os.remove(f)
                    raise Exception("Failed to download "+f+". You may retry \
                            to run the script to complete the download.")
        existing_files.add(f)


parser = argparse.ArgumentParser(description = """Downloads the databases necessary
        for the annotation pipeline""")

parser.add_argument("--download_ko", nargs="?", const="databases/ko/", default=False,
        help="download ko file definition, necessary to do this before --load-kegg")
parser.add_argument("--download_refseq", nargs="?", const="databases/refseq/", default=False,
        help="download ko file definition, necessary to do this before --load-kegg")


args = vars(parser.parse_args())

threads = []


if args.get("download_refseq", False):
    download_dir = args["download_refseq"]
    thread = threading.Thread(target=download_refseq, args=[download_dir])
    threads.append(thread)
    thread.start()


if args.get("download_ko", False):
    print("Download ko")
    pass


# wait until completion
for thread in threads:
    thread.join()
