#!/usr/bin/env python

# get complete genome sequences from ncbi using eutils
# TODO replace prints by sys.stdout/err.write
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 01.2015
# ---------------------------------------------------------------------------

from Bio import Entrez, SeqIO


Entrez.email = "trestan.pillonel@unil.ch"

"""
plsmids: not working?
handle_plasmids = Entrez.elink(dbfrom="genome", db="nuccore", id=one_genome_id, term="srcdb+ddbj/embl/genbank[prop] AND gene+in+plasmid[prop]")
record_plasmids = Entrez.read(handle_plasmids)


if len(record_plasmids[0]["LinkSetDb"][0]["Link"]) == 0:
    print "No plasmid seq for %s" % one_genome_id
else:
    linked_plasmids = [link["Id"] for link in record_plasmids[0]["LinkSetDb"][0]["Link"]]
    print "Plasmid(s):", linked_plasmids
    genome_record_id_list += linked_plasmids
"""

def download_wgs(ncbi_taxon_id):
    links = multiple_wgs_links(ncbi_taxon_id)
    for link in links:
        download_one_wgs(link)

def download_one_wgs(wgs_link):
        handle = Entrez.elink(dbfrom="nuccore", db="nuccore", id=wgs_link)
        record = Entrez.read(handle)
        # get the list of link
        try:
            sublinks = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
        except IndexError:
            print "No link for %s" % wgs_link
            return None
        #print "sublinks", sublinks

#        if len(sublinks) > 500:
#            print "More than 500 contigs, aborting download for id %s" % wgs_link
#            return None

        # for each sublink (contig, scaffold,...), get record and append the sequence to a single file
        output_handle = open("%s.gbk" % wgs_link, "a")
        no_sequences = False

        
        for seq_link in sublinks:
            if no_sequences:
                print "No sequences fo record %s" % wgs_link
                output_handle.close()
                
                import os
                os.remove("%s.gbk" % wgs_link)
                
                break

            while not handle:
                print i
                if i == 10:
                    print 'reached max iteration number, %s could not be downloaded' % record_id
                    return
                try:
                    handle = Entrez.efetch(db="nucleotide", id=seq_link, rettype="gb", retmode="text")
                except (urllib2.URLError, urllib2.HTTPError) as e:
                    print 'url error, trying again...'
                    time.sleep(1)
                    i+=1                
            
            seq_records = list(SeqIO.parse(handle, "genbank"))
            # in case multiple record for a single link (shouldn't append???)
            for record in seq_records:
                print record.name
                print record.description
                if record.seq.count("N") == len(record.seq):
                    no_sequences = True
                    break
                else:
                    SeqIO.write(record, output_handle, "genbank")
        output_handle.close()    


def get_wgs_links(one_species_link):

        # get all WGS linked to this species
        handle = Entrez.elink(dbfrom="genome", db="nuccore", id=one_species_link, term="wgs[prop]")
        record = Entrez.read(handle)
        if len(record[0]["LinkSetDb"][0]["Link"]) == 0:
            print "No WGS genome seq for %s" % one_species_link
            return False
        else:
            linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
            #print "WGS genome(s):", linked
            return linked


def multiple_wgs_links(ncbi_taxon):
    #handle = Entrez.esearch(db="genome", term="klebsiella+pneumoniae[orgn]")
    #txid570[Organism:exp]
    handle = Entrez.esearch(db="genome", term="txid%s[Organism:exp]" % ncbi_taxon)
    record = Entrez.read(handle)

    # get genome overview id
    genome_id_list = record["IdList"]

    # elink.fcgi?dbfrom=genome&db=nuccore&id=1076

    # get whole genomes only
    genome_record_id_list = []
    for one_genome_id in genome_id_list:
        one_genome_ids = get_wgs_links(one_genome_id)
        if one_genome_ids:
            genome_record_id_list += one_genome_ids
        
    return genome_record_id_list
        
        
def get_complete_genomes_data(ncbi_taxon):
    import time
    import eutils
    import urllib2
    #handle = Entrez.esearch(db="genome", term="klebsiella+pneumoniae[orgn]")
    #txid570[Organism:exp]
    handle = Entrez.esearch(db="genome", term="txid%s[Organism:exp]" % ncbi_taxon)
    record = Entrez.read(handle)
    print record
    # get genome overview id
    genome_id_list = record["IdList"]

    # elink.fcgi?dbfrom=genome&db=nuccore&id=1076

    # get whole genomes only
    genome_record_id_list = []
    for one_genome_id in genome_id_list:
        print "considering genome ID", one_genome_id
        # srcdb+ddbj/embl/genbank[prop] AND 
        handle = Entrez.elink(dbfrom="genome", db="nuccore", id=one_genome_id, term="gene+in+genomic[prop] OR gene+in+chromosome[prop]")


        record = Entrez.read(handle)
        print "all links", record[0]["LinkSetDb"][0]["Link"]
        
        if len(record[0]["LinkSetDb"][0]["Link"]) == 0:
            print "No whole genome seq for %s" % one_genome_id

        else:
            print "Whole genome data for %s " % one_genome_id
            print record
            
            linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
            print "Complete genome(s):", linked
            genome_record_id_list += linked

    n = 1
    for record_id in genome_record_id_list:
        print "ID:", record_id, "(%s out of %s)" % (n, len(genome_record_id_list))
        n+=1
        # get assembly
        i = 0
        handle_assembly = None
        while not handle_assembly:
            print i
            if i == 10:
                print 'reached max iteration number, %s could not be downloaded' % record_id
                return
            try:
                handle_assembly = Entrez.elink(dbfrom="nuccore", db="assembly", id=record_id)
            except (urllib2.URLError, urllib2.HTTPError) as e:
                print 'url error, trying again...'
                time.sleep(1)
                i+=1                
                
        record_assembly = Entrez.read(handle_assembly)
        try:
            assembly_link = record_assembly[0]["LinkSetDb"][0]["Link"][0]["Id"]
            print "assembly link", assembly_link
            # assembly 2 genome + plasmids
            handle_sequences = Entrez.elink(dbfrom="assembly", db="nuccore", id=assembly_link, term="srcdb+ddbj/embl/genbank[prop]")
            record_sequences =  Entrez.read(handle_sequences)
            wgs = False
            for i in range(0, len(record_sequences[0]["LinkSetDb"])):
                if record_sequences[0]["LinkSetDb"][i]['LinkName'] == 'assembly_nuccore_wgsmaster':
                    wgs = True
            if wgs:
                print "wgs link"
                continue
            
            sequences_links = [link["Id"] for link in record_sequences[0]["LinkSetDb"][0]["Link"]]            
        except:
            print "No assembly link, searching bioproject for %s..." % record_id

            # first: check if not WGS:
            handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="gb", retmode="text")
            seq_records = list(SeqIO.parse(handle, "genbank"))


            print "annot", seq_records[0].annotations
            if "wgs" in seq_records[0].annotations or "WGS" in seq_records[0].annotations:
                print "Wgs, not downloading %s" % record_id, seq_records[0].description
                continue                
            elif isinstance(seq_records[0].annotations, dict) and len(seq_records[0].annotations.values()) > 1:
                tag = False
                for i in seq_records[0].annotations.values():
                    print "annot", i
                    if isinstance(i, list):
                        if "wgs" in i or "WGS" in i:
                            print "Wgs, not downloading %s" % record_id, seq_records[0].description
                            tag = True
                    elif isinstance(i, str):
                        if "wgs" in i or "WGS" in i:
                            print "Wgs, not downloading %s" % record_id, seq_records[0].description
                            tag = True
                    elif isinstance(i, int):
                        continue
                    else:
                        print "problem with:", i
                if tag:
                    continue

            handle_bioproject = Entrez.elink(dbfrom="nuccore", db="bioproject", id=record_id)
            record_bioproject = Entrez.read(handle_bioproject)

            # get bioproject link
            bioproject_link = record_bioproject[0]["LinkSetDb"][0]["Link"][0]["Id"]

            # as bioprojects are not linked with genbank, get link(s) to refseq (for both plasmids and chromosomes)
            handle_refseq = Entrez.elink(dbfrom="bioproject", db="nuccore", id=bioproject_link)
            record_refseq =  Entrez.read(handle_refseq)

            refseq_links = [link["Id"] for link in record_refseq[0]["LinkSetDb"][0]["Link"]]

            print "refseqlinks", refseq_links
            # for each refseq id, get link to genbank entry (which contain full record)
            sequences_links = []
            for link in refseq_links:
                handle_genbank = Entrez.elink(dbfrom="nuccore", db="nuccore", id=link, term="srcdb+ddbj/embl/genbank[prop]")
                record_genbank = Entrez.read(handle_genbank)
                sequences_links.append(record_genbank[0]["LinkSetDb"][1]["Link"][0]["Id"])

        print "Sequences links:", sequences_links

        if len(sequences_links) == 0:
            print "PROBLEM WITH ASSEMBLY"
            print record_sequences
        #import time
        #time.sleep(4)
        for one_id in sequences_links:
            eutils.get_genomic_data(one_id)





if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--taxon_id',type=str,help="get wgs link from taxonomic id")
    parser.add_argument("-w",'--wgs',type=str,help="download one wgs link")
    parser.add_argument("-c",'--complete_genomes',type=str,help="taxonomic id (get complete genomes only")

    parser.add_argument("-a",'--all_genomes', type=str,help="taxonomic id (get all genomes (wgs + complete)")



    
    args = parser.parse_args()
    if args.all_genomes:
        download_wgs(args.all_genomes)
        get_complete_genomes_data(args.all_genomes)

        import sys
        sys.exit()



    if args.taxon_id:
        print "getting link to wgs records for taxon %s ..." % args.taxon_id
        all_links = multiple_wgs_links(args.taxon_id)
        for i in all_links:
            print i

    if args.wgs:
        print "downloading record %s..." % args.wgs
        download_one_wgs(args.wgs)
   
    if args.complete_genomes:
        print "downloading complete genomes from taxon %s..." % args.complete_genomes
        get_complete_genomes_data(args.complete_genomes)
