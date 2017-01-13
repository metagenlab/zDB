#!/usr/bin/python

def orthogroup2uniprot_hits(biodb, orthogroup, max_n_hits=10, exclude_phylum="Chlamydiae",mysql_host="localhost",
                            mysql_user="root", mysql_pwd="estrella3", mysql_db="blastnr"):

    # get max n hits for each homolog of <orthogroup> from  other phyla
    # remove redundancy
    # return fasta file
    import MySQLdb
    conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                user=mysql_user, # your username
                                passwd=mysql_pwd, # your password
                                db=mysql_db) # name of the data base
    cursor = conn.cursor()

    sql = 'select t2.locus_tag,subject_accession from biosqldb.orthology_detail_%s t1 ' \
          ' inner join custom_tables.locus2seqfeature_id_%s t2 ' \
          ' on t1.locus_tag=t2.locus_tag ' \
          ' inner join blast_swissprot_%s t3 on t2.seqfeature_id=t3.seqfeature_id ' \
          ' inner join blastnr_taxonomy as t4 on t3.subject_taxid=t4.taxon_id ' \
          ' where t1.orthogroup="%s" and t4.phylum!="%s" order by locus_tag,hit_number;' % (biodb,
                                                                                            biodb,
                                                                                            biodb,
                                                                                            orthogroup,
                                                                                            exclude_phylum)
    cursor.execute(sql,)
    data = cursor.fetchall()

    locus2top_n_hits = {}

    for row in data:
        if row[0] not in locus2top_n_hits:
            locus2top_n_hits[row[0]] = [row[1]]
        else:
            if len(locus2top_n_hits[row[0]]) < max_n_hits:
                locus2top_n_hits[row[0]].append(row[1])

    return locus2top_n_hits

def swissprot_accession2fasta(accession_list):

    import urllib2
    from urllib2 import URLError
    from Bio import SeqIO
    import StringIO

    accession_list_req = ','.join(accession_list)

    link = "http://www.uniprot.org/uniprot/?query=%s&format=fasta" % (accession_list_req)
    print link
    req = urllib2.Request(link)

    page = urllib2.urlopen(req)

    handle = StringIO.StringIO(page.read().decode('utf-8'))

    print 'download ok'

    records = [i for i in SeqIO.parse(handle, 'fasta')]

    return records

def orthogroup2locus_and_sequences(biodb,orthogroup):

    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select locus_tag, translation from orthology_detail_%s where orthogroup="%s"' % (biodb, orthogroup)

    data = server.adaptor.execute_and_fetchall(sql,)

    return data


def orthogroup2alignment_closest(orthogroup, biodb, max_n_hits_uniprot=2, exclude_phylum="chlamydiae", outname="out.fa"):

    locus2uniprot_accession_list = orthogroup2uniprot_hits(biodb, orthogroup, max_n_hits=max_n_hits_uniprot, exclude_phylum=exclude_phylum)

    uniprot_accession_list = []
    for locus in locus2uniprot_accession_list:
        uniprot_accession_list+=locus2uniprot_accession_list[locus]

    print 'n hits:', len(uniprot_accession_list)
    print uniprot_accession_list


    uniprot_sequence_records = swissprot_accession2fasta(uniprot_accession_list)

    ortho_sequences = orthogroup2locus_and_sequences(biodb, orthogroup)

    with open(outname, 'w') as f:
        for locus in ortho_sequences:
            f.write(">%s\n%s\n" % (locus[0], locus[1]))
        for record in uniprot_sequence_records:
            name = record.name.split('|')[1]
            f.write(">%s\n%s\n" % (name, str(record.seq)))


def aafasta2phylogeny(aa_fasta):
    import os
    import shell_command

    align_name = aa_fasta.split('.')[0] + '_mafft.fa'

    cmd_mafft = 'mafft --anysymbol --amino --auto --maxiterate 1000 %s > %s' % (aa_fasta, align_name)

    out, err, code = shell_command.shell_command(cmd_mafft)

    if code != 0:
        raise(err)

    output_prefix = aa_fasta.split('.')[0]
    output_tree_name = os.path.join('RAxML_result.%s' % output_prefix)
    output_shtree_name = os.path.join('shtest_%s' % output_prefix)
    cmd_raxml = 'raxml -m PROTGAMMALG -p 12345 -s %s -n %s -c 4 -T 8;' \
          'raxml -f J -m PROTGAMMALG -s %s -p 12345 -t %s -n %s -T 8' % (align_name,
                                                                     output_prefix,
                                                                     align_name,
                                                                     output_tree_name,
                                                                     output_shtree_name)

    out, err, code = shell_command.shell_command(cmd_raxml)
    print out, err, code

orthogroup2alignment_closest("group_444", "chlamydia_04_16", 2, "Chlamydiae", outname="group_444_swiss_homologs.faa")

aafasta2phylogeny("group_444_swiss_homologs.faa")