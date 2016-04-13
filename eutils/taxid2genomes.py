#!/usr/bin/env python

# get complete genome sequences from any ncbi taxon
# 1. identify assemblies related to the taxon ID
# 2. get the assembly ftp link from the assemblies metedata
# 3. download the assemblies from the ftp
# NOTE: download RefSeq assembly first, if no RefSeq assembly, download GenBank assembly
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 01.2015
# ---------------------------------------------------------------------------

from Bio import Entrez, SeqIO
import os
import re
from ftplib import FTP
import download_from_ftp

Entrez.email = "trestan.pillonel@unil.ch"

def get_ncbi_genome(path, destination):
    ftp=FTP('ftp.ncbi.nih.gov')
    ftp.login("anonymous","trestan.pillonel@unil.ch")
    download_from_ftp.download_whole_directory(ftp, path, destination, recursive=False)

def get_complete_genomes_data(ncbi_taxon, complete = False):

    handle = Entrez.esearch(db="assembly", term="txid%s[Organism:exp]" % ncbi_taxon, retmax=100000)
    record = Entrez.read(handle)
    f = open('ncbi_genome_download.log', 'w')
    header = 'AssemblyName\tLastMajorReleaseAccession\tTaxid\tSpeciesName\tOrganism\t' \
             'AssemblyStatus\tNCBIReleaseDate\tSubmitterOrganization\tftpPath\tPartialGenomeRepresentation\tCoverage\tRefSeq\tGenbank\tSimilarity\tAnomalousAssembly\n'
    f.write(header)
    local_dir = os.getcwd()

    for i, one_assembly in enumerate(record['IdList']):
        print 'assembly %s (%s/%s)' % (one_assembly, i, len(record['IdList']))
        handle_assembly = Entrez.esummary(db="assembly",id=one_assembly)
        assembly_record = Entrez.read(handle_assembly, validate=False)


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
        if complete:
            if AssemblyStatus not in ['Complete Genome']: # , 'Chromosome'
                print 'status:', AssemblyStatus
                continue

        NCBIReleaseDate = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['NCBIReleaseDate']
        SpeciesName = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['SpeciesName']
        SubmitterOrganization = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['SubmitterOrganization']
        AssemblyName = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']
        RefSeq = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['RefSeq']
        Genbank = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']
        Similarity = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Similarity']
        PartialGenomeRepresentation = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PartialGenomeRepresentation']
        Coverage = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Coverage']

        out_dir = local_dir+'/%s' % LastMajorReleaseAccession
        try:
            os.mkdir(out_dir)
        except:
            continue
        # get ftp link in meta data string
        print assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta']
        print assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']

        if 'anomalous' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            anomalous = True
        else:
            anomalous = False

        # check if RefSeq assembly was supressed
        if 'suppressed_refseq' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            print 'Supressed RefSeq record, downloading GenBank record'
            ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
            get_ncbi_genome(ftp_path, out_dir)
        elif 'latest' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            print 'Downloading RefSeq record'
            try:
                ftp_path = re.findall('<FtpPath type="RefSeq">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
                ftp_path = '/' + ftp_path
                get_ncbi_genome(ftp_path, out_dir)
            except IndexError:
                print 'RefSeq link could not be identified, downloading GenBank record'
                ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
                get_ncbi_genome(ftp_path, out_dir)

        # check if the assembly was replaced
        elif 'replaced' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:

            i = int(LastMajorReleaseAccession.split('.')[1]) + 1
            success = False
            # identifiy new version of the assembly by increasing version number
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
        line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (AssemblyName,
                                                                             LastMajorReleaseAccession,
                                                                             Taxid,
                                                                             SpeciesName,
                                                                             Organism,
                                                                             AssemblyStatus,
                                                                             NCBIReleaseDate,
                                                                             SubmitterOrganization,
                                                                             ftp_path,
                                                                             PartialGenomeRepresentation,
                                                                             Coverage,
                                                                             RefSeq,
                                                                             Genbank,
                                                                             Similarity,
                                                                             anomalous)

        f.write(line.encode("utf-8"))
    f.close()


def get_sample_data(accession):

    handle1 = Entrez.esearch(db="nuccore", term=accession)
    record1 = Entrez.read(handle1)

    ncbi_id = record1['IdList'][0]

    handle2 = Entrez.elink(dbfrom="nuccore", db="biosample", id=ncbi_id)
    record2 = Entrez.read(handle2)

    id = record2[0]['LinkSetDb'][0]['Link'][0]['Id']

    handle3 = Entrez.esummary(db='biosample',id=id, retmode='xml')
    record = Entrez.read(handle3, validate=False)
    for i in record['DocumentSummarySet']['DocumentSummary'][0]:
        print "%s\t%s" % (i, record['DocumentSummarySet']['DocumentSummary'][0][i])


def accession2assembly(accession):

    handle1 = Entrez.esearch(db="nuccore", term=accession)
    record1 = Entrez.read(handle1)

    print record1

    ncbi_id = record1['IdList'][-1]

    print "ncbi_id", ncbi_id

    handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=ncbi_id)
    record = Entrez.read(handle)

    print record

    f = open('ncbi_genome_download.log', 'w')
    header = 'AssemblyName\tLastMajorReleaseAccession\tTaxid\tSpeciesName\tOrganism\t' \
             'AssemblyStatus\tNCBIReleaseDate\tSubmitterOrganization\tftpPath\tPartialGenomeRepresentation\tCoverage\tRefSeq\tGenbank\tSimilarity\n'
    f.write(header)
    local_dir = os.getcwd()

    handle_assembly = Entrez.esummary(db="assembly",id=record[0]['LinkSetDb'][0]['Link'][0]['Id'])
    assembly_record = Entrez.read(handle_assembly)

    LastMajorReleaseAccession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['LastMajorReleaseAccession']
    Taxid = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Taxid']
    Organism = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
    AssemblyStatus = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus']
    NCBIReleaseDate = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['NCBIReleaseDate']
    SpeciesName = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['SpeciesName']
    SubmitterOrganization = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['SubmitterOrganization']
    AssemblyName = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']
    RefSeq = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['RefSeq']
    Genbank = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']
    Similarity = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Similarity']
    PartialGenomeRepresentation = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PartialGenomeRepresentation']
    Coverage = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Coverage']

    out_dir = local_dir+'/%s' % LastMajorReleaseAccession
    os.mkdir(out_dir)

    # get ftp link in meta data string
    if 'latest' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
        ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]

        get_ncbi_genome(ftp_path, out_dir)

    elif 'replaced' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:

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
    line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (AssemblyName,
                                                                         LastMajorReleaseAccession,
                                                                         Taxid,
                                                                         SpeciesName,
                                                                         Organism,
                                                                         AssemblyStatus,
                                                                         NCBIReleaseDate,
                                                                         SubmitterOrganization,
                                                                         ftp_path,
                                                                         PartialGenomeRepresentation,
                                                                         Coverage,
                                                                         RefSeq,
                                                                         Genbank,
                                                                         Similarity)

    f.write(line.encode("utf-8"))



def bioproject2assemblies(accession, complete=False):

    handle1 = Entrez.esearch(db="bioproject", term=accession)
    record1 = Entrez.read(handle1)

    ncbi_id = record1['IdList'][0]

    handle = Entrez.elink(dbfrom="bioproject", db="assembly", id=ncbi_id)
    record = Entrez.read(handle)

    id_list = [i['Id'] for i in record[0]['LinkSetDb'][0]['Link']]
    print id_list


    local_dir = os.getcwd()
    f = open('ncbi_genome_download.log', 'w')
    header = 'AssemblyName\tLastMajorReleaseAccession\tTaxid\tSpeciesName\tOrganism\t' \
             'AssemblyStatus\tNCBIReleaseDate\tSubmitterOrganization\tftpPath\tPartialGenomeRepresentation\tCoverage\tRefSeq\tGenbank\tSimilarity\tAnomalousAssembly\n'
    f.write(header)

    for i, one_assembly in enumerate(id_list):
        print 'assembly %s (%s/%s)' % (one_assembly, i, len(id_list))
        handle_assembly = Entrez.esummary(db="assembly",id=one_assembly)
        assembly_record = Entrez.read(handle_assembly, validate=False)

        LastMajorReleaseAccession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['LastMajorReleaseAccession']

        Taxid = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Taxid']
        Organism = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
        AssemblyStatus = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus']
        if complete:
            if AssemblyStatus not in ['Complete Genome', 'Chromosome']:
                print 'status:', AssemblyStatus
                continue

        NCBIReleaseDate = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['NCBIReleaseDate']
        SpeciesName = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['SpeciesName']
        SubmitterOrganization = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['SubmitterOrganization']
        AssemblyName = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']
        RefSeq = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['RefSeq']
        Genbank = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']
        Similarity = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Similarity']
        PartialGenomeRepresentation = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PartialGenomeRepresentation']
        Coverage = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Coverage']

        out_dir = local_dir+'/%s' % LastMajorReleaseAccession
        os.mkdir(out_dir)

        # get ftp link in meta data string
        print assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta']
        print assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']

        if 'anomalous' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            anomalous = True
        else:
            anomalous = False

        # check if RefSeq assembly was supressed
        if 'suppressed_refseq' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            print 'Supressed RefSeq record, downloading GenBank record'
            ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
            get_ncbi_genome(ftp_path, out_dir)
        elif 'latest' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            print 'Downloading RefSeq record'
            try:
                ftp_path = re.findall('<FtpPath type="RefSeq">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
                ftp_path = '/' + ftp_path
                get_ncbi_genome(ftp_path, out_dir)
            except IndexError:
                print 'RefSeq link could not be identified, downloading GenBank record'
                ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
                get_ncbi_genome(ftp_path, out_dir)

        # check if the assembly was replaced
        elif 'replaced' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:

            i = int(LastMajorReleaseAccession.split('.')[1]) + 1
            success = False
            # identifiy new version of the assembly by increasing version number
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
        line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (AssemblyName,
                                                                             LastMajorReleaseAccession,
                                                                             Taxid,
                                                                             SpeciesName,
                                                                             Organism,
                                                                             AssemblyStatus,
                                                                             NCBIReleaseDate,
                                                                             SubmitterOrganization,
                                                                             ftp_path,
                                                                             PartialGenomeRepresentation,
                                                                             Coverage,
                                                                             RefSeq,
                                                                             Genbank,
                                                                             Similarity,
                                                                             anomalous)

        f.write(line.encode("utf-8"))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--taxon_id',type=str,help="get wgs link from taxonomic id")
    parser.add_argument("-a",'--accession',type=str,help="get assembly of a specific genbank/refseq accession")
    parser.add_argument("-c", '--complete', action="store_true", help="complete genome only (default = False)")
    parser.add_argument("-s",'--sample',type=str,help="get sample data from accession id")
    parser.add_argument("-b",'--bioproject',type=str,help="get all assemblies from bioproject")

    args = parser.parse_args()
    if args.taxon_id:
        get_complete_genomes_data(ncbi_taxon= args.taxon_id, complete = args.complete)
    if args.accession:
        accession2assembly(accession = args.accession)
    if args.sample:
        get_sample_data(accession=args.sample)
    if args.bioproject:
        bioproject2assemblies(accession=args.bioproject)
