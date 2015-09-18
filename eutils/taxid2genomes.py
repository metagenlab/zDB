#!/usr/bin/env python

# get complete genome sequences from ncbi using eutils
# TODO replace prints by sys.stdout/err.write
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 01.2015
# ---------------------------------------------------------------------------

from Bio import Entrez, SeqIO
import urllib2
import time
import os
Entrez.email = "trestan.pillonel@unil.ch"

import ftplib
import re
from ftplib import FTP
import download_from_ftp


def get_ncbi_genome(path, destination):
    ftp=FTP('ftp.ncbi.nih.gov')
    ftp.login("anonymous","trestan.pillonel@unil.ch")
    download_from_ftp.download_whole_directory(ftp, path, destination, recursive=False)


def get_complete_genomes_data(ncbi_taxon):

    handle = Entrez.esearch(db="assembly", term="txid%s[Organism:exp]" % ncbi_taxon, retmax=100000)
    record = Entrez.read(handle)
    f = open('ncbi_genome_download.log', 'w')
    header = 'AssemblyName\tLastMajorReleaseAccession\tTaxid\tSpeciesName\tOrganism\t' \
             'AssemblyStatus\tNCBIReleaseDate\tSubmitterOrganization\tftpPath\n'
    f.write(header)
    local_dir = os.getcwd()
    for one_assembly in record['IdList']:

        handle_assembly = Entrez.esummary(db="assembly",id=one_assembly)
        assembly_record = Entrez.read(handle_assembly)


        '''
        # check if it is the latest available version of the genome, if not, download the latest
        if 'latest' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            pass
        else:
            handle = Entrez.esearch(db="assembly",
                                    term=assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank'],
                                    retmax=100000)
            record = Entrez.read(handle)
            print record['IdList'][0]
            handle_assembly = Entrez.esummary(db="assembly",id=record['IdList'][0])
            assembly_record = Entrez.read(handle_assembly)
        '''
        LastMajorReleaseAccession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['LastMajorReleaseAccession']

        Taxid = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Taxid']
        Organism = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
        AssemblyStatus = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus']
        NCBIReleaseDate = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['NCBIReleaseDate']
        SpeciesName = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['SpeciesName']
        SubmitterOrganization = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['SubmitterOrganization']
        AssemblyName = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']


        out_dir = local_dir+'/%s' % LastMajorReleaseAccession
        os.mkdir(out_dir)

        # get ftp link in meta data string
        if 'latest' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]

            get_ncbi_genome(ftp_path, out_dir)


        elif 'replaced' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            #print 'Assembly %s was replaced: not downloading' % LastMajorReleaseAccession
            #ftp_path = '/genomes/all/%s_%s' % (AssemblyName)
            #print assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym'], '######################'


            i = int(LastMajorReleaseAccession.split('.')[1]) + 1
            success = False
            while not success:
                if i < 5:
                    try:
                        ftp_path = '/genomes/all/%s_%s' % (LastMajorReleaseAccession.split('.')[0] + '.%s' % i,
                                                           AssemblyName)
                        print ftp_path
                        get_ncbi_genome(ftp_path, out_dir)

                        success = True
                    except:
                        i+=i
                else:
                    print 'could not download assembly %s' % LastMajorReleaseAccession
                    break
        else:
            'Unkown status for %s, not downloading' % LastMajorReleaseAccession
        line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (AssemblyName, LastMajorReleaseAccession, Taxid,
                                                       SpeciesName, Organism, AssemblyStatus, NCBIReleaseDate,
                                                       SubmitterOrganization, ftp_path)
        f.write(line)


        #ftp_folder = LastMajorReleaseAccession + '_' + AssemblyName
        #print ftp_folder





    f.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--taxon_id',type=str,help="get wgs link from taxonomic id")



    args = parser.parse_args()
    get_complete_genomes_data(args.taxon_id)