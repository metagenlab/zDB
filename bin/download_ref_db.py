#!/usr/bin/env python
import argparse
import ftplib
import os
import re
import time


def download_refseq(download_dir, n_retry=10):
    ftp = ftplib.FTP("ftp.ncbi.nih.gov")
    ftp.login("anonymous")
    ftp.cwd("/refseq/release/complete")

    # check if there are any pre-existing files in the download directory
    # to avoid downloading them a second time
    # NOTE: this assumes that the files are either completely downloaded
    # or not present (i.e. removed from the disk in case of problem)
    existing_files = set()
    for current_file in os.listdir(download_dir):
        if os.path.isfile(download_dir + "/" + current_file):
            existing_files.add(current_file)

    # it would probably be better to use startsWith and check the last
    # letters to make sure it ends with faa.gz, instead of using regular
    # expression
    nr_re = re.compile(r"complete.nonredundant_protein.(\d)*.protein.faa.gz")
    nr_filelist = [
        i
        for i in ftp.nlst()
        if re.match(nr_re, i) is not None and i not in existing_files
    ]

    for f in nr_filelist:
        failed = 0
        complete = False
        while not complete:
            output_file = open(download_dir + "/" + f, "wb")
            try:
                ftp.retrbinary("RETR " + f, output_file.write)
                output_file.close()
                complete = True
                time.sleep(5)
            except Exception:
                output_file.close()
                time.sleep(60)
                failed += 1
                if failed == n_retry:
                    os.remove(f)
                    raise Exception(
                        "Failed to download "
                        + f
                        + ". You may retry \
                            to run the script to complete the download."
                    )
        existing_files.add(f)


parser = argparse.ArgumentParser(
    description="""Downloads the databases necessary
        for the annotation pipeline"""
)

parser.add_argument(
    "--download_ko",
    nargs="?",
    const="databases/ko/",
    default=False,
    help="download ko file definition, necessary to do this before --load-kegg",
)
parser.add_argument(
    "--download_refseq",
    nargs="?",
    const="databases/refseq/",
    default=False,
    help="download ko file definition, necessary to do this before --load-kegg",
)


args = vars(parser.parse_args())


if args.get("download_refseq", False):
    download_dir = args["download_refseq"]
    download_refseq(download_dir)


if args.get("download_ko", False):
    print("Download ko")
    pass
