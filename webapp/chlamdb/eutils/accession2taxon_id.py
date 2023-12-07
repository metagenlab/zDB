#!/usr/bin/env python

from Bio import Entrez
import urllib.request

Entrez.email = "trestan.pillonel@unil.ch"


def old_taxon_id2new_taxon_id(old_taxon_id):
    # handle = Entrez.esearch(db="taxonomy", term="%s[uid]" % old_taxon_id)
    handle = Entrez.esummary(db="taxonomy", id=old_taxon_id)
    record1 = Entrez.read(handle)
    print('record1', record1)
    print(record1[0]['AkaTaxId'])
    try:
        if record1[0]['Status'] == 'merged':
            ncbi_id = record1[0]['AkaTaxId']
            return ncbi_id
        else:
            raise ('Unkown status %s' % record1[0]['Status'])
    except:
        return False


def gi(ncbi_term, database="nuccore", retmax=20):

    handle = Entrez.esearch(db=database, term=ncbi_term, retmax=retmax)

    record = Entrez.read(handle)

    return record["IdList"]


def species_name2taxon_id(species_name):
    from socket import error as SocketError
    import errno
    handle = Entrez.esearch(db="taxonomy", term=species_name)
    record1 = Entrez.read(handle)
    try:
        ncbi_id = record1['IdList'][0]
        return ncbi_id
    except:
        return False


def gi2taxon_id(gi_list, database="protein"):
    from socket import error as SocketError
    import errno

    print('babababa')
    try:
        handle = Entrez.esummary(
            db=database, id=','.join(gi_list), retmax=len(gi_list))
    except SocketError as e:
        if e.errno != errno.ECONNRESET:
            print('error connexion with %s' % ','.join(gi_list))
        else:
            import time
            print('connexion error, trying again...')
            time.sleep(60)
            gi2taxon_id(gi_list, database=database)

    gi2taxon = {}
    print('ok')
    print(isinstance(gi_list, list))
    print(len(gi_list))
    if isinstance(gi_list, list) and len(gi_list) == 1:
        record = Entrez.read(handle, validate=False)
        print(record)
        gi2taxon[record['AccessionVersion']] = record['TaxId']
        return gi2taxon
    else:
        record = Entrez.parse(handle, validate=False)
        try:

            for i in record:
                print('gg', i)
                # i['Gi']
                # print i
                gi2taxon[i['AccessionVersion']] = i['TaxId']
        except RuntimeError:
            gi2taxon_id(gi_list, database=database)

        return gi2taxon


def gi2protein_accession(gi, database="nuccore"):
    from socket import error as SocketError
    import errno

    if database != 'protein':
        try:
            handle = Entrez.elink(dbfrom=database, id=gi, db='protein')
        except SocketError as e:
            if e.errno != errno.ECONNRESET:
                print('error connexion with %s' % ','.join(gi))
            else:
                import time
                print('connexion error, trying again...')
                time.sleep(60)
                gi2protein_accession(gi, database=database)
        try:
            record = Entrez.read(handle, validate=False)
            print('record', record)
            try:
                gi_prot = record[0]['LinkSetDb'][0]['Link'][0]['Id']
            except IndexError:
                return None
            handle = Entrez.esummary(db='protein', id=gi_prot)
            record = Entrez.read(handle, validate=False)

        except RuntimeError:
            gi2protein_accession(gi, database=database)
    else:
        handle = Entrez.esummary(db='protein', id=gi)
        record = Entrez.read(handle, validate=False)
    return record[0]['AccessionVersion']


def accession2gbk(accession, db='protein'):
    from Bio import SeqIO

    try:
        handle = Entrez.esearch(db=db, term="%s" % (accession), retmax=1)
        record = Entrez.read(handle)
    except urllib.request.HTTPError:
        import time
        time.sleep(5)
        return accession2gbk(accession, db=db)
    try:
        gi_id = record['IdList'][0]
    except KeyError:
        print('try again!')
        return locus_tag2protein_accession(locus_tag)
    if db == 'gene':
        handle3 = Entrez.efetch(db="gene", id=gi_id, rettype="xml")
        record3 = Entrez.read(handle3, validate=False)
        gi_id = record3[0]['Entrezgene_xtra-iq'][1]['Xtra-Terms_value']

    try:
        handlexx = Entrez.efetch(db=db,
                                 id="%s" % gi_id,
                                 rettype="gb", retmode="text")

        record_gbk = SeqIO.read(handlexx, 'genbank')
        return record_gbk
    except:
        return None


def locus_tag2protein_accession(locus_tag, db='protein'):
    from socket import error as SocketError
    import errno
    # print 'locus', locus_tag
    handle = Entrez.esearch(db=db, term="%s" % (locus_tag), retmax=1)
    record = Entrez.read(handle)
    print(record)
    try:
        gi_id = record['IdList'][0]
    except KeyError:
        print('try again!')
        return locus_tag2protein_accession(locus_tag)
    if db == 'gene':
        handle3 = Entrez.efetch(db="gene", id=gi_id, rettype="xml")
        record3 = Entrez.read(handle3, validate=False)
        gi_id = record3[0]['Entrezgene_xtra-iq'][1]['Xtra-Terms_value']

    handle2 = Entrez.esummary(db='protein', id="%s" % gi_id)

    try:
        record = Entrez.read(handle2, validate=False)
        print(record)

        if record[0]['Status'] == 'dead':
            if 'ReplacedBy' in record[0]:
                return [record[0]['ReplacedBy'], record[0]['TaxId']]
            else:
                return None
        elif record[0]['Status'] == 'suppressed':
            return [record[0]['AccessionVersion'], record[0]['TaxId']]

        else:
            return [record[0]['AccessionVersion'], record[0]['TaxId']]
    except RuntimeError:
        print('try again!')
        return locus_tag2protein_accession(locus_tag)


def accession2taxon_id(ncbi_id_list, db="protein"):
    '''
    :param genbank/refseq protein accession
    :return: dictionnary protein_id2taxon_id
    '''
    # print ncbi_id_list
    # print
    # print 'n unique taxons', len(set(ncbi_id_list))

    gi_list = gi(','.join(ncbi_id_list), db, retmax=len(ncbi_id_list))

    # print 'n unique gi', len(set(gi_list))

    # get link from genbank 2 refseq
    handle = Entrez.elink(dbfrom=db, db="taxonomy", id=','.join(gi_list))
    record = Entrez.read(handle)

    # print 'len record', len(record)

    try:
        # if no result, return None
        # print 'record length:', len(record)
        refseq_ids = [i['LinkSetDb'][0]['Link'][0]['Id'] for i in record]
        # print "taxon_ids", refseq_ids
        # print 'n taxons:', len(ncbi_id_list), 'n gi:', len(gi_list), 'n taxon ids', len(refseq_ids)
        if len(ncbi_id_list) != len(refseq_ids):
            print('multiple matches for %s' % ncbi_id_list)
            print(ncbi_id_list)
            print(gi_list)
            print(refseq_ids)
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

    # print "gi_list", ','.join(ncbi_gi_list)
    # get link from genbank 2 refseq
    match = False
    while not match:
        # print ncbi_gi_list
        handle = Entrez.efetch(db=db, id=','.join(
            ncbi_gi_list), rettype="gb", retmode="text")

        try:
            records = [i for i in SeqIO.parse(handle, "genbank")]
            match = True
        except IncompleteRead:
            print('problem trying to fetch data for list:')
            print(','.join(ncbi_gi_list))

    res_list = []

    # print len(records)
    for record in records:
        # print record
        # print dir(record)
        tmp = {}
        # accession2description[record.id]
        # print 'record', record.id
        # print dir(record)
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

    print("ncbi_id_list", ncbi_id_list)
    gi_list = gi(','.join(ncbi_id_list), db)
    print("gi_list", ','.join(gi_list))
    # get link from genbank 2 refseq
    handle = Entrez.efetch(db=db, id=gi_list, rettype="gb", retmode="text")

    records = [i for i in SeqIO.parse(handle, "genbank")]

    res_list = []

    # print 'records', len(records)
    for record in records:
        tmp = {}
        # accession2description[record.id]
        # print 'record', record.id
        # print dir(record)
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
    print('baba')
    accession2taxon = accession2taxon_id(accession_list, database)
    taxon_id_list = accession2taxon.values()

    # print 'taxon_id list',

    classif = sequence_id2scientific_classification.taxon_id2scientific_classification(
        taxon_id_list)

    return classif, accession2taxon


if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--seq_id_genbank', default=False,
                        type=str, help="genbank2refseq", nargs="+")
    parser.add_argument("-p", '--full_path', action="store_true",
                        help="get full classification path")
    parser.add_argument("-e", '--description',
                        help="get description", action="store_true")
    parser.add_argument("-g", '--gi2protein_accesson',
                        help="gi2protein acession", action="store_true")
    parser.add_argument("-d", '--ncbi_database', default="nucleotide",
                        type=str, help="database to search (protein/nucleotide/...)")

    args = parser.parse_args()

    # species_name2taxon_id("Gemmata obscuriglobus baba")
    # import sys
    # sys.exit()

    if args.full_path:
        print('full path!')
        taxon2path, accession2taxon = accession2full_taxonomic_path(
            args.seq_id_genbank, args.ncbi_database)
        print('ok')
        # {'superkingdom': 'Bacteria',
        # 'no rank': 'Chlamydia/Chlamydophila group',
        # 'phylum': 'Chlamydiae',
        # 'class': 'Chlamydiia'
        # 'order': 'Chlamydiales',
        # 'family': 'Chlamydiaceae',
        # 'species': 'Chlamydia muridarum',
        # 'genus': 'Chlamydia',
        # }
        # print "superkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies"
        if taxon2path:
            for accession in taxon2path:
                data = taxon2path[accession]
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
                print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (accession2taxon.keys()[accession2taxon.values().index(accession)],
                                                              accession,
                                                              superkingdom,
                                                              phylum,
                                                              class_,
                                                              order,
                                                              family,
                                                              genus,
                                                              species))
        else:
            print(args.seq_id_genbank)
    elif args.description:
        # print 'description!'
        accession2description(args.seq_id_genbank, args.ncbi_database)

    else:

        if args.gi2protein_accesson:
            print('ok', gi2protein_accession(
                args.seq_id_genbank, args.ncbi_database))
        else:
            # print 'normal!'
            dico = accession2taxon_id(args.seq_id_genbank, args.ncbi_database)
            for accession in dico:
                print('%s\t%s' % (accession, dico[accession]))
