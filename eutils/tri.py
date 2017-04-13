def get_complete_genomes_data_old(ncbi_taxon):

    import eutils

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
        handle = Entrez.elink(dbfrom="genome", db="nuccore", id=one_genome_id, term="gene+in+chromosome[prop]")
        record = Entrez.read(handle)

        print "all links", record[0]["LinkSetDb"][0]["Link"]

        if len(record[0]["LinkSetDb"][0]["Link"]) == 0:
            print "No whole genome seq for %s" % one_genome_id
        else:
            linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
            print "Complete genome(s):", linked
            genome_record_id_list += linked
    print "Complete list:", genome_record_id_list

    for record_id in genome_record_id_list:
        print "ID:", record_id
        eutils.get_genomic_data(record_id)

    # get wgs only wgs[prop]
    #print "WGS"
    #handle = Entrez.elink(dbfrom="genome", db="nuccore", id=one_genome_id, term="wgs[prop]")
    #record = Entrez.read(handle)
    #print record

    #handle = Entrez.elink(dbfrom="nucleotide", db="nuccore", id=ncbi_accession, term="contig")
    #record = Entrez.read(handle)
    #print record[0]
    #print(record[0]["LinkSetDb"][0]["LinkName"])

    #print "all linked WGS sequences:", record[0]["LinkSetDb"][0]["Link"]

    #linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]

    #print linked
    #print len(linked)

    #for link in linked:
    #    print "link", link
    #    handle = Entrez.efetch(db="nucleotide", id=link, rettype="gb", retmode="text")
    #    seq_records = list(SeqIO.parse(handle, "genbank"))
    #    for record in seq_records:
    #        print record.name
    #        print record.description

    """
    handle = Entrez.efetch(db="nucleotide", id=ncbi_accession, rettype="gbwithparts", retmode="text")
    seq_records = list(SeqIO.parse(handle, "genbank"))


    seq_records = list(SeqIO.parse(handle, "genbank"))
    for record in seq_records:
        print "writing record %s..." % record.name
        SeqIO.write(record, "%s.gbk" % record.name, "genbank")
        SeqIO.write(record, "%s.fna" % record.name, "fasta")
        gbk2faa(record, "%s.faa" % record.name)
        gbk2ffn(record, "%s.ffn" % record.name)
    handle.close()
    """


def get_complete_genomes_data2(ncbi_taxon):

    import eutils

    #handle = Entrez.esearch(db="genome", term="klebsiella+pneumoniae[orgn]")
    #txid570[Organism:exp]
    handle = Entrez.search(db="genome", term="txid%s[Organism:exp]" % ncbi_taxon)
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
        handle = Entrez.elink(dbfrom="genome", db="nuccore", id=one_genome_id, term="gene+in+chromosome[prop]")


        record = Entrez.read(handle)
        print "all links", record[0]["LinkSetDb"][0]["Link"]

        if len(record[0]["LinkSetDb"][0]["Link"]) == 0:
            print "No whole genome seq for %s" % one_genome_id
        else:
            print "Whole genome data for %s " % one_genome_id
            print record

            """
            handle_plasmids = Entrez.elink(dbfrom="genome", db="nuccore", id=one_genome_id, term="srcdb+ddbj/embl/genbank[prop] AND gene+in+plasmid[prop]")
            record_plasmids = Entrez.read(handle_plasmids)


            if len(record_plasmids[0]["LinkSetDb"][0]["Link"]) == 0:
                print "No plasmid seq for %s" % one_genome_id
            else:
                linked_plasmids = [link["Id"] for link in record_plasmids[0]["LinkSetDb"][0]["Link"]]
                print "Plasmid(s):", linked_plasmids
                genome_record_id_list += linked_plasmids
            """

            linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
            print "Complete genome(s):", linked
            genome_record_id_list += linked



    print "Complete list:", genome_record_id_list, len(genome_record_id_list)

    for record_id in genome_record_id_list:
        print "ID:", record_id

        # get assembly
        handle_assembly = Entrez.elink(dbfrom="nuccore", db="assembly", id=record_id)
        record_assembly = Entrez.read(handle_assembly)
        assembly_link = record_assembly[0]["LinkSetDb"][0]["Link"][0]["Id"]
        print "assembly link", assembly_link

        # assembly 2 genome + plasmids
        handle_sequences = Entrez.elink(dbfrom="assembly", db="nuccore", id=assembly_link, term="srcdb+ddbj/embl/genbank[prop]")
        record_sequences =  Entrez.read(handle_sequences)

        sequences_links = [link["Id"] for link in record_sequences[0]["LinkSetDb"][0]["Link"]]

        print "Sequences links:", sequences_links

        if len(sequences_links) == 0:
            print "PROBLEM WITH ASSEMBLY"
            print record_sequences
        #import time
        #time.sleep(4)
        for one_id in sequences_links:
            eutils.get_genomic_data(one_id)

