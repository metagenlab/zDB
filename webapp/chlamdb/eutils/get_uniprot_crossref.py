#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"

def gi(ncbi_term, database="nuccore", retmax=20):

    handle = Entrez.esearch(db=database, term=ncbi_term, retmax=retmax)
    record = Entrez.read(handle)

    return record["IdList"]


def gi2taxon_id(gi_list, database="protein"):
    from socket import error as SocketError
    import errno

    try:
        handle = Entrez.esummary(db=database, id=','.join(gi_list), retmax=len(gi_list))
    except SocketError as e:
        if e.errno != errno.ECONNRESET:
            print ('error connexion with %s' % ','.join(gi_list))
        else:
            import time
            print ('connexion error, trying again...')
            time.sleep(60)
            gi2taxon_id(gi_list, database=database)
    record = Entrez.parse(handle, validate=False)
    gi2taxon = {}
    try:
        for i in record:
            gi2taxon[i['Gi']] = i['TaxId']
    except RuntimeError:
        gi2taxon_id(gi_list, database=database)

    return gi2taxon



def accession2taxon_id(ncbi_id_list, db="protein"):
    '''
    :param genbank/refseq protein accession
    :return: dictionnary protein_id2taxon_id
    '''
    #print ncbi_id_list
    #print
    #print 'n unique taxons', len(set(ncbi_id_list))

    gi_list = gi(','.join(ncbi_id_list), db, retmax=len(ncbi_id_list))

    #print 'n unique gi', len(set(gi_list))

    # get link from genbank 2 refseq
    handle = Entrez.elink(dbfrom=db, db="taxonomy", id=','.join(gi_list))
    record = Entrez.read(handle)


    try:
        # if no result, return None
        #print 'record length:', len(record)
        refseq_ids = [i['LinkSetDb'][0]['Link'][0]['Id'] for i in record]
        #print "taxon_ids", refseq_ids
        #print 'n taxons:', len(ncbi_id_list), 'n gi:', len(gi_list), 'n taxon ids', len(refseq_ids)
        if len(ncbi_id_list) != len(refseq_ids):
            print (ncbi_id_list)
            print (gi_list)
            print (refseq_ids)
            import sys
            sys.exit()
        dico = dict(zip(ncbi_id_list, refseq_ids))
        return dico
    except IndexError:
        return None

def gi2description(ncbi_gi_list, db="protein"):
    '''
    :param genbank/refseq protein accession
    :return: dictionnary protein_id2description
    '''
    from Bio import SeqIO
    import re

    #print "gi_list", ','.join(ncbi_gi_list)
    # get link from genbank 2 refseq
    match = False
    while not match:
        handle = Entrez.efetch(db=db, id=','.join(ncbi_gi_list), rettype="gb", retmode="text")

        try:
            records = [i for i in SeqIO.parse(handle, "genbank")]
            match = True
        except IncompleteRead:
            print ('problem trying to fetch data for list:')
            print (','.join(ncbi_gi_list))
        
    res_list = []

    for record in records:
        #print record
        #print dir(record)
        tmp = {}
        #accession2description[record.id]
        #print 'record', record.id
        #print dir(record)
        tmp['description'] = re.sub('\"', '', record.description).split('[')[0]
        tmp['taxonomy'] = record.annotations["taxonomy"]
        tmp['source'] = record.annotations['source']
        tmp['division'] = record.annotations['data_file_division']
        res_list.append(tmp)
    return dict(zip(ncbi_gi_list, res_list))

def accession2description(ncbi_id_list, db="protein"):
    '''
    :param genbank/refseq protein accession
    :return: dictionnary protein_id2description
    '''
    from Bio import SeqIO
    import re

    print ("ncbi_id_list", ncbi_id_list)
    gi_list = gi(','.join(ncbi_id_list), db)
    print ("gi_list", ','.join(gi_list))
    # get link from genbank 2 refseq
    handle = Entrez.efetch(db=db, id=gi_list, rettype="gb", retmode="text")

    records = [i for i in SeqIO.parse(handle, "genbank")]

    res_list = []

    #print len(records)
    for record in records:
        tmp = {}
        #accession2description[record.id]
        #print 'record', record.id
        #print dir(record)
        tmp['description'] = re.sub('\"', '', record.description).split('[')[0]
        tmp['taxonomy'] = record.annotations["taxonomy"]
        tmp['source'] = record.annotations['source']
        tmp['division'] = record.annotations['data_file_division']
        res_list.append(tmp)
    return dict(zip(ncbi_id_list, res_list))

    '''
    try:
        # if no result, return None
        refseq_ids = [ i['LinkSetDb'][0]['Link'][0]['Id'] for i in record]
        print "refseq_ids", refseq_ids
        dico = dict(zip(ncbi_id_list, refseq_ids))
        return dico
    except IndexError:
        return None
    '''


def accession2full_taxonomic_path(accession_list, database="nucleotide"):
    import sequence_id2scientific_classification

    taxon_id_list = accession2taxon_id(accession_list, database).values()

    #print 'taxon_id list',

    classif = sequence_id2scientific_classification.taxon_id2scientific_classification(taxon_id_list)

    return classif

if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--seq_id_genbank', default=False, type=str, help="genbank2refseq", nargs="+")
    parser.add_argument("-p",'--full_path', action="store_true", help="get full classification path")
    parser.add_argument("-e", '--description', help="get description", action="store_true")
    parser.add_argument("-d",'--ncbi_database', default="nucleotide", type=str, help="database to search (protein/nucleotide/...)")

    args = parser.parse_args()


    if args.full_path:
        dico = accession2full_taxonomic_path(args.seq_id_genbank, args.ncbi_database)
        # {'superkingdom': 'Bacteria',
        # 'no rank': 'Chlamydia/Chlamydophila group',
        # 'phylum': 'Chlamydiae',
        # 'class': 'Chlamydiia'
        # 'order': 'Chlamydiales',
        # 'family': 'Chlamydiaceae',
        # 'species': 'Chlamydia muridarum',
        # 'genus': 'Chlamydia',
        # }
        #print "superkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies"
        if dico:
            for accession in dico:
                data = dico[accession]
                try:
                    superkingdom = data["superkingdom"]
                except:
                    superkingdom = "-"
                try:
                    phylum = data["phylum"]
                except:
                    phylum = "-"
                try:
                    class_ = data["class"]
                except:
                    class_ = "-"
                try:
                    order = data["order"]
                except:
                    order = "-"
                try:
                    family = data["family"]
                except:
                    family = "-"
                try:
                    genus = data["genus"]
                except:
                    genus = "-"
                try:
                    species = data["species"]
                except:
                    species = "-"

                print ("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (superkingdom, phylum, class_, order, family, genus, species))
        else:
            print (args.seq_id_genbank)
    elif args.description:
        accession2description(args.seq_id_genbank, args.ncbi_database)

    else:
        print (accession2taxon_id(args.seq_id_genbank, args.ncbi_database))
