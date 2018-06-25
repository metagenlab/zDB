#!/usr/bin/env python

def process_tag(tag):
    return tag.split('}')[-1]

def get_UPI(seq):
    for element in seq:
        if element.tag == '{http://model.picr.ebi.ac.uk}UPI':
            return element.text

def get_hit_attributes(hit):
    accession = ''
    version = ''
    taxon_id = ''
    db_name = ''
    for element in hit:
        if element.tag == '{http://model.picr.ebi.ac.uk}accession':
            accession = element.text
        if element.tag == '{http://model.picr.ebi.ac.uk}accessionVersion':
            version = element.text
        if element.tag == '{http://model.picr.ebi.ac.uk}databaseName':
            db_name = element.text
        if element.tag == '{http://model.picr.ebi.ac.uk}taxonId':
            taxon_id = element.text
    return {"%s.%s" % (accession, version) : [db_name, taxon_id]}

def accession2exact_matches(sequence, target_databases):

    '''
    Givent an input AA sequence and target(s) database name(s), return:
     - the uniparc accession of the sequence (if exists)
     - a dictionary with accession(s) of identical sequence(s) and their taxon ID and source database.
       (Accession.version keys)

    Return None if no identical squence was found.

    :param sequence: input AA sequence
    :param target_databases: Input database name (see http://www.ebi.ac.uk/Tools/picr/)

    '''


    import urllib2
    import xml.etree.cElementTree as ElementTree

    database_string = '&database=' .join(target_databases)

    link = "http://www.ebi.ac.uk/Tools/picr/rest/getUPIForSequence?sequence=%s&database=%s&includeattributes=true" % (sequence,
                                                                                                                      database_string)
    print link

    req = urllib2.Request(link)
    try:
        page = urllib2.urlopen(req)
        tree = ElementTree.parse(page)
    except:
        import time
        print 'connexion problem, trying again...'
        time.sleep(60)

    db2seq = {}
    root = tree.getroot()

    seq = root.find('{http://www.ebi.ac.uk/picr/AccessionMappingService}getUPIForSequenceReturn')
    if seq is None:
        return None
    UPI = get_UPI(seq)
    identical_seqs = seq.findall('{http://model.picr.ebi.ac.uk}identicalCrossReferences')
    for seq in identical_seqs:
        db2seq.update(get_hit_attributes(seq))
    return UPI, db2seq


def fasta_corresp(fasta_file, target_database, n_keep=1):
    from Bio import SeqIO
    import sys

    print 'keep', n_keep

    with open(fasta_file, 'r') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:

            picr = accession2exact_matches(record.seq,
                                           target_database)

            if picr is None:
                sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % (record.name, 'None', 'None', 'None', record.seq))
            else:
                uniparc_accession, matches = picr
                database2count = {}
                for accession in matches:
                    if matches[accession][0] not in database2count:
                        database2count[matches[accession][0]] = 1
                    else:
                        if database2count[matches[accession][0]] < n_keep:
                            database2count[matches[accession][0]] += 1
                        else:

                            break
                    sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % (record.name,
                                                               uniparc_accession,
                                                               accession,
                                                               matches[accession][1],
                                                               record.seq))
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", '--protein_seq', type=str, help="Protein sequence")
    parser.add_argument("-d", '--database', type=str, help="Target database(s): 'REFSEQ', 'TREMBL', ...", nargs='+', default= ['TREMBL', 'SWISSPROT'])
    parser.add_argument("-f", '--fasta_file', type=str, help="Fasta file")
    parser.add_argument("-k", '--keep', type=int, help="Number of hit(s) to keep (default: 1)", default=1)

    args = parser.parse_args()

    if args.protein_seq and args.fasta_file:
        raise(IOError('Input either a fasta file or a protein seqience, not both!'))
    elif args.protein_seq:
        picr = accession2exact_matches(args.protein_seq,
                                       args.database)
        if picr is not None:
            uniparc_accession, matches = picr
            print uniparc_accession, matches
    else:
        if len(args.database) > 1:
            raise(IOError('Fasta file match is only possible for a single database!'))
        else:
            fasta_corresp(args.fasta_file, args.database, n_keep=args.keep)
