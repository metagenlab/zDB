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
    
    return download_from_ftp.download_whole_directory(ftp, path, destination, recursive=False)


def get_taxi2assembly_accession(ncbi_taxon,
                              complete = False,
                              representative =False,
                              refseq_only=False,
                              genbank=False):
    s_term = 'txid%s[Organism:exp] (AND all[filter] NOT anomalous[filter]) '
    if complete:
        s_term+=' AND "complete genome"[filter]'
    if representative:
        s_term+=' AND "representative genome"[filter]'
    if refseq_only:
        s_term+=' AND "latest refseq"[filter]'
    else:
        s_term+= ' AND latest[filter]'
    handle = Entrez.esearch(db="assembly", term=s_term % ncbi_taxon, retmax=100000)
    record = Entrez.read(handle)
    for i, one_assembly in enumerate(record['IdList']):
        try:
            handle_assembly = Entrez.esummary(db="assembly", id=one_assembly)
        except urllib2.URLError:
            import time
            urlok = False
            while not urlok:
                print ('connexion problem, waiting 10s...')
                time.sleep(10)
                try:
                    handle_assembly = Entrez.esummary(db="assembly",id=one_assembly)
                    urlok = True
                except:
                    urlok = False
        try:
            assembly_record = Entrez.read(handle_assembly, validate=False)
        except Bio.Entrez.Parser.CorruptedXMLError:
            print ('assembly:', one_assembly)
        #print assembly_record['DocumentSummarySet']
        AssemblyName = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']
        RefSeq = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['RefSeq']
        Genbank = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']
        try:
            NCBIReleaseDate = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['NCBIReleaseDate']
        except KeyError:
            NCBIReleaseDate = '-'
        contig_count = re.findall('<Stat category="contig_count" sequence_tag="all">([^<]*)<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0]
        status = re.findall('<assembly-status>([^<]*)<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0]
        repres = re.findall('representative-status>([^<]*)<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0]
        length = re.findall('Stat category="total_length" sequence_tag="all">([^<]*)<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0]
        #print assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta']
        print ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (ncbi_taxon,
                                          AssemblyName,
                                          RefSeq,
                                          Genbank,
                                          contig_count,
                                          status,
                                          repres,
                                              length))


def get_complete_genomes_data(ncbi_taxon,
                              complete = False,
                              representative =False,
                              refseq_only=False,
                              genbank=False):

    s_term = 'txid%s[Organism:exp] (AND all[filter] NOT anomalous[filter]) '
    if complete:
        s_term+=' AND "complete genome"[filter]'
    if representative:
        s_term+=' AND "representative genome"[filter]'
    if refseq_only:
        s_term+=' AND "latest refseq"[filter]'
    else:
        s_term+= ' AND latest[filter]'
    handle = Entrez.esearch(db="assembly", term=s_term % ncbi_taxon, retmax=100000)

    record = Entrez.read(handle)
    f = open('ncbi_genome_download.log', 'w')
    nn = open("assembly_log.tx", "w")
    header = 'AssemblyName\tLastMajorReleaseAccession\tTaxid\tSpeciesName\tOrganism\t' \
             'AssemblyStatus\tNCBIReleaseDate\tSubmitterOrganization\tftpPath\tPartialGenomeRepresentation\tCoverage\tRefSeq\tGenbank\tSimilarity\tAnomalousAssembly\n'
    f.write(header)
    local_dir = os.getcwd()

    for i, one_assembly in enumerate(record['IdList']):
        print ("%s/%s" % (i, len(record['IdList'])))
        try:
            handle_assembly = Entrez.esummary(db="assembly", id=one_assembly)
        except urllib2.URLError:
            import time
            urlok = False
            while not urlok:
                print ('connexion problem, waiting 10s...')
                time.sleep(10)
                try:
                    handle_assembly = Entrez.esummary(db="assembly",id=one_assembly)
                    urlok = True
                except:
                    urlok = False
        try:
            assembly_record = Entrez.read(handle_assembly, validate=False)
        except Bio.Entrez.Parser.CorruptedXMLError:
            print ('assembly:', one_assembly)

        LastMajorReleaseAccession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['LastMajorReleaseAccession']
        SubDate = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['SubmissionDate']
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
        PartialGenomeRepresentation = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PartialGenomeRepresentation']
        Coverage = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Coverage']
        try:
            Biosource = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Biosource']['InfraspeciesList'][0]['Sub_value']
        except:
            try:
                Biosource = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Biosource']['Isolate']
            except:
                print (assembly_record['DocumentSummarySet']['DocumentSummary'])
                import sys
                sys.exit()

        nn.write('assembly %s (%s/%s), %s, %s, %s, %s, %s, %s, %s\n' % (one_assembly, i, len(record['IdList']),
                                                                        LastMajorReleaseAccession,
                                                                        AssemblyStatus,
                                                                        Organism,
                                                                        SpeciesName,
                                                                        Biosource,
                                                                        SubDate,
                                                                        assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']))

        out_dir = local_dir+'/%s' % LastMajorReleaseAccession
        try:
            os.mkdir(out_dir)
        except:
            continue
        # get ftp link in meta data string
        #print assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta']
        #print assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']

        # check if RefSeq assembly was supressed ==> if yes get the genbank link
        if 'suppressed_refseq' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            print ('Supressed RefSeq record, downloading GenBank record')
            ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
            get_ncbi_genome(ftp_path, out_dir)
        elif 'latest' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            print ('Downloading RefSeq record')
            try:
                ftp_path = re.findall('<FtpPath type="RefSeq">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
                ftp_path = '/' + ftp_path
                status = get_ncbi_genome(ftp_path, out_dir)

                if status == False:
                    if genbank:
                        print ('RefSeq link is empty, downloading GenBank record')
                        ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
                        get_ncbi_genome(ftp_path, out_dir)
                    else:
                        print ('no RefSeq record for assembly %s' % (one_assembly))
                        continue
            except IndexError:

                if genbank:
                    print ('RefSeq link could not be identified, downloading GenBank record')
                    ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
                    get_ncbi_genome(ftp_path, out_dir)
                else:
                    print ('no RefSeq record for assembly %s' % (one_assembly))
                    continue
                try:
                    print ('RefSeq ftp link could not be found, downloading GenBank record')
                    #print assembly_record['DocumentSummarySet']['DocumentSummary'][0]
                    ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
                    get_ncbi_genome(ftp_path, out_dir)
                except:
                    print ('no ftp link found for assembly id:', one_assembly)
                    os.rmdir(out_dir)
                    continue

        # check if the assembly was replaced
        elif 'replaced' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            print ('replaced version-search for latest-------------------------')
            i = int(LastMajorReleaseAccession.split('.')[1]) + 1
            success = False
            # identifiy latest version of the assembly by increasing version number
            while not success:
                if i < 5:
                    try:
                        ftp_path = '/genomes/all/%s_%s' % (LastMajorReleaseAccession.split('.')[0] + '.%s' % i,
                                                           AssemblyName)
                        print (ftp_path)
                        get_ncbi_genome(ftp_path, out_dir)

                        success = True
                    except:
                        i+=i
                else:
                    print ('could not download assembly %s' % LastMajorReleaseAccession)
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
                                                                             Similarity
                                                                             )

        f.write(line)
        os.chdir(local_dir)

    f.close()
    nn.close()


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
        print ("%s\t%s" % (i, record['DocumentSummarySet']['DocumentSummary'][0][i]))


def assembly_accession2assembly(accession):
    handle1 = Entrez.esearch(db="assembly", term=accession)
    record1 = Entrez.read(handle1)


    ncbi_id = record1['IdList'][-1]

    local_dir = os.getcwd()
    try:
        handle_assembly = Entrez.esummary(db="assembly", id=ncbi_id)
    except:

        print ('link to assembly could not be found, trying to get it strating from taxon_id...')

        handle_assembly = Entrez.esummary(db="nucleotide",id=ncbi_id)
        record1 = Entrez.read(handle_assembly)
        taxon_id = record1[0]['TaxId']
        term = "txid%s[Organism:noexp]" % taxon_id
        handle1 = Entrez.esearch(db="assembly", term=term)
        record1 = Entrez.read(handle1)
        assembly_id = record1['IdList'][0]
        handle_assembly = Entrez.esummary(db="assembly",id=assembly_id)


    assembly_record = Entrez.read(handle_assembly, validate=False)

    LastMajorReleaseAccession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['LastMajorReleaseAccession']
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
    PartialGenomeRepresentation = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PartialGenomeRepresentation']
    Coverage = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Coverage']

    out_dir = local_dir+'/%s' % LastMajorReleaseAccession
    os.mkdir(out_dir)

    # get ftp link in meta data string
    if 'latest' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
        print (assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])
        ftp_path = re.findall('<FtpPath type="RefSeq">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]

        get_ncbi_genome(ftp_path, out_dir)

    elif 'replaced' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:

        i = int(LastMajorReleaseAccession.split('.')[1]) + 1
        success = False
        while not success:
            if i < 5:
                try:
                    ftp_path = '/genomes/all/%s_%s' % (LastMajorReleaseAccession.split('.')[0] + '.%s' % i,
                                                       AssemblyName)
                    print (ftp_path)
                    get_ncbi_genome(ftp_path, out_dir)

                    success = True
                except:
                    i+=i
            else:
                print ('could not download assembly %s' % LastMajorReleaseAccession)
                break
    else:
        'Unkown status for %s, not downloading' % LastMajorReleaseAccession


def accession2assembly(accession):

    handle1 = Entrez.esearch(db="nuccore", term=accession)
    record1 = Entrez.read(handle1)

    print ('noccore hit:', record1)

    ncbi_id = record1['IdList'][-1]

    print ("ncbi_id", ncbi_id)

    handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=ncbi_id)
    record = Entrez.read(handle)

    print ('links:', record)

    f = open('ncbi_genome_download.log', 'w')
    header = 'AssemblyName\tLastMajorReleaseAccession\tTaxid\tSpeciesName\tOrganism\t' \
             'AssemblyStatus\tNCBIReleaseDate\tSubmitterOrganization\tftpPath\tPartialGenomeRepresentation\tCoverage\tRefSeq\tGenbank\tSimilarity\n'
    f.write(header)
    local_dir = os.getcwd()
    try:
        handle_assembly = Entrez.esummary(db="assembly",id=record[0]['LinkSetDb'][0]['Link'][0]['Id'])
    except:

        print ('link to assembly could not be found, trying to get it strating from taxon_id...')

        handle_assembly = Entrez.esummary(db="nucleotide",id=ncbi_id)
        record1 = Entrez.read(handle_assembly)
        taxon_id = record1[0]['TaxId']
        term = "txid%s[Organism:noexp]" % taxon_id
        handle1 = Entrez.esearch(db="assembly", term=term)
        record1 = Entrez.read(handle1)
        assembly_id = record1['IdList'][0]
        handle_assembly = Entrez.esummary(db="assembly",id=assembly_id)


    assembly_record = Entrez.read(handle_assembly, validate=False)

    LastMajorReleaseAccession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['LastMajorReleaseAccession']
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
    PartialGenomeRepresentation = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PartialGenomeRepresentation']
    Coverage = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Coverage']

    out_dir = local_dir+'/%s' % LastMajorReleaseAccession
    os.mkdir(out_dir)

    # get ftp link in meta data string
    if 'latest' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
        print (assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])
        ftp_path = re.findall('<FtpPath type="RefSeq">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]

        get_ncbi_genome(ftp_path, out_dir)

    elif 'replaced' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:

        i = int(LastMajorReleaseAccession.split('.')[1]) + 1
        success = False
        while not success:
            if i < 5:
                try:
                    ftp_path = '/genomes/all/%s_%s' % (LastMajorReleaseAccession.split('.')[0] + '.%s' % i,
                                                       AssemblyName)
                    print (ftp_path)
                    get_ncbi_genome(ftp_path, out_dir)

                    success = True
                except:
                    i+=i
            else:
                print ('could not download assembly %s' % LastMajorReleaseAccession)
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


def download_assembly(assembly_gi, out_dir, complete=False):

        handle_assembly = Entrez.esummary(db="assembly",id=assembly_gi)
        assembly_record = Entrez.read(handle_assembly, validate=False)

        LastMajorReleaseAccession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['LastMajorReleaseAccession']

        Taxid = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Taxid']
        Organism = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
        AssemblyStatus = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus']
        if complete:
            if AssemblyStatus not in ['Complete Genome', 'Chromosome']:
                print ('status:', AssemblyStatus)
                return False
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
        PartialGenomeRepresentation = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PartialGenomeRepresentation']
        Coverage = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Coverage']

        out_dir = out_dir+'/%s' % LastMajorReleaseAccession
        os.mkdir(out_dir)

        # get ftp link in meta data string
        print (assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])
        print (assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList'])

        if 'anomalous' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            anomalous = True
        else:
            anomalous = False

        # check if RefSeq assembly was supressed
        if 'suppressed_refseq' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            print ('Supressed RefSeq record, downloading GenBank record')
            ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
            get_ncbi_genome(ftp_path, out_dir)
        elif 'latest' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]['PropertyList']:
            print ('Downloading RefSeq record')
            try:
                ftp_path = re.findall('<FtpPath type="RefSeq">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
                print ('\nftp path: %s\n' % ftp_path)
                ftp_path = '/' + ftp_path
                get_ncbi_genome(ftp_path, out_dir)
            except IndexError:
                print ('RefSeq link could not be identified, downloading GenBank record')
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
                        print (ftp_path)
                        get_ncbi_genome(ftp_path, out_dir)

                        success = True
                    except:
                        i+=i
                else:
                    print ('could not download assembly %s' % LastMajorReleaseAccession)
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

        return line

def bioproject2assemblies(accession, complete=False):

    handle1 = Entrez.esearch(db="bioproject", term=accession)
    record1 = Entrez.read(handle1)

    ncbi_id = record1['IdList'][0]

    handle = Entrez.elink(dbfrom="bioproject", db="assembly", id=ncbi_id)
    record = Entrez.read(handle)

    id_list = [i['Id'] for i in record[0]['LinkSetDb'][0]['Link']]
    print (id_list)


    local_dir = os.getcwd()
    f = open('ncbi_genome_download.log', 'w')
    header = 'AssemblyName\tLastMajorReleaseAccession\tTaxid\tSpeciesName\tOrganism\t' \
             'AssemblyStatus\tNCBIReleaseDate\tSubmitterOrganization\tftpPath\tPartialGenomeRepresentation\tCoverage\tRefSeq\tGenbank\tSimilarity\tAnomalousAssembly\n'
    f.write(header)

    for i, one_assembly in enumerate(id_list):
        print ('assembly %s (%s/%s)' % (one_assembly, i, len(id_list)))
        dw = download_assembly(one_assembly, local_dir, complete)
        if dw:
            f.write(dw.encode("utf-8"))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--taxon_id',type=str,help="get wgs link from taxonomic id")
    parser.add_argument("-a",'--accession',type=str,help="get assembly of a specific genbank/refseq accession")
    parser.add_argument("-c", '--complete', action="store_true", help="complete genome only (default = False)")
    parser.add_argument("-s",'--sample',type=str,help="get sample data from accession id")
    parser.add_argument("-b",'--bioproject',type=str,help="get all assemblies from bioproject")
    parser.add_argument("-r",'--representative',action="store_true", help="download only representative genomes")
    parser.add_argument("-g", '--genbank', action="store_true", help="Download only Genbank record if the assembly is not in RefSeq")
    parser.add_argument("-aa", '--assembly_accession', type=str, help="Download based on assembly accession")
    parser.add_argument("-ta", '--taxid2assembly_accession', action="store_true", help="get correspondance table from taxid2accession")

    args = parser.parse_args()
    if args.taxon_id:
        if args.taxid2assembly_accession:
            get_taxi2assembly_accession(ncbi_taxon=args.taxon_id,
                                        complete=args.complete,
                                        representative=args.representative,
                                        genbank=args.genbank)
        else:
            get_complete_genomes_data(ncbi_taxon= args.taxon_id,
                                  complete = args.complete,
                                  representative = args.representative,
                                  genbank=args.genbank)
    if args.accession:
        accession2assembly(accession=args.accession)
    if args.sample:
        get_sample_data(accession=args.sample)
    if args.bioproject:
        bioproject2assemblies(accession=args.bioproject)
    if args.assembly_accession:
        assembly_accession2assembly(args.assembly_accession)

