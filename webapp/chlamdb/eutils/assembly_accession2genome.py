#!/usr/bin/env python

# download assembly from genbank or refseq based on accession number
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 10.2018
# ---------------------------------------------------------------------------

from Bio import Entrez, SeqIO
import os
import re
from ftplib import FTP
import download_from_ftp
import taxid2genomes

Entrez.email = "trestan.pillonel@unil.ch"


def assembly_accession2assembly(accession,
                                only_gbk=False,
                                only_fna=True):

    handle1 = Entrez.esearch(db="assembly", term=accession)
    record1 = Entrez.read(handle1)

    ncbi_id = record1['IdList'][-1]

    local_dir = os.getcwd()
    try:
        handle_assembly = Entrez.esummary(db="assembly", id=ncbi_id)
    except Exception:

        print('link to assembly could not be found, trying to get it strating from taxon_id...')

        handle_assembly = Entrez.esummary(db="nucleotide", id=ncbi_id)
        record1 = Entrez.read(handle_assembly)
        taxon_id = record1[0]['TaxId']
        term = "txid%s[Organism:noexp]" % taxon_id
        handle1 = Entrez.esearch(db="assembly", term=term)
        record1 = Entrez.read(handle1)
        assembly_id = record1['IdList'][0]
        handle_assembly = Entrez.esummary(db="assembly", id=assembly_id)

    assembly_record = Entrez.read(handle_assembly, validate=False)

    LastMajorReleaseAccession = assembly_record['DocumentSummarySet'][
        'DocumentSummary'][0]['LastMajorReleaseAccession']
    Taxid = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Taxid']
    Organism = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
    AssemblyStatus = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus']
    try:
        NCBIReleaseDate = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['NCBIReleaseDate']
    except KeyError:
        NCBIReleaseDate = '-'
    SpeciesName = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['SpeciesName']
    SubmitterOrganization = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['SubmitterOrganization']
    AssemblyName = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']
    RefSeq = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['RefSeq']
    Genbank = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']
    Similarity = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Similarity']
    PartialGenomeRepresentation = assembly_record['DocumentSummarySet'][
        'DocumentSummary'][0]['PartialGenomeRepresentation']
    Coverage = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Coverage']

    out_dir = local_dir + '/%s' % LastMajorReleaseAccession
    os.mkdir(out_dir)

    # get ftp link in meta data string
    if 'latest' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
        print(assembly_record['DocumentSummarySet']
              ['DocumentSummary'][0]['Meta'])
        try:
            ftp_path = re.findall(
                '<FtpPath type="RefSeq">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
        except IndexError:
            print('RefSeq ftp link could not be found, downloading GenBank record')
            # print assembly_record['DocumentSummarySet']['DocumentSummary'][0]
            ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<',
                                  assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
        taxid2genomes.get_ncbi_genome(ftp_path, out_dir,
                                      only_gbk=only_gbk,
                                      only_fna=only_fna)

    elif 'replaced' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:

        i = int(LastMajorReleaseAccession.split('.')[1]) + 1
        success = False
        while not success:
            if i < 5:
                try:
                    ftp_path = '/genomes/all/%s_%s' % (LastMajorReleaseAccession.split('.')[0] + '.%s' % i,
                                                       AssemblyName)
                    print(ftp_path)
                    get_ncbi_genome(ftp_path, out_dir)

                    success = True
                except Exception:
                    i += i
            else:
                print('could not download assembly %s' %
                      LastMajorReleaseAccession)
                break
    else:
        'Unkown status for %s, not downloading' % LastMajorReleaseAccession


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--assembly_accession',
                        type=str, help="NCBI assembly accession")
    parser.add_argument("-f", '--only_fna', action='store_true',
                        help="Download only fna file")
    parser.add_argument("-g", '--only_gbk', action='store_true',
                        help="Download only gbff file")

    args = parser.parse_args()
    assembly_accession2assembly(args.assembly_accession,
                                only_gbk=args.only_gbk,
                                only_fna=args.only_fna)
